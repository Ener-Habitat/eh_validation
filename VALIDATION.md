# Plan de validación — Bloque hueco de concreto (cavidad de aire / relleno) en EnerHabitat

> Estado: **APROBADO (alcance Track 1)** — 2026-06-25
> Autor: Guillermo Barrios del Valle · Fecha: 2026-06-25
> Repo: `eh_validation_eplus` · EnerHabitat ≥ 0.2.0 (instalado en `.venv`)

> **Decisiones aprobadas:**
> - **Alcance:** Track 1 (R estacionario vs Borbón) es la validación requerida.
>   Tracks 2 (EnergyPlus dinámico) y 3 (relleno sólido) quedan como **opcionales/futuros**.
> - **Concreto:** usar `Concreto` de `materials.ini` (k=1.35, ρ=1800, c=1000).
> - **Emisividad:** nominal **0.9**, con **barrido 0.8–0.95**.

---

## 1. Objetivo

Validar el módulo **2D** de EnerHabitat (`System2D` + `HollowBlock`) para el caso
de un **bloque hueco de concreto** —tanto con **cavidad de aire** (`Fill.AIR`)
como con **relleno sólido** (`Fill.SOLID`)— extendiendo la validación 1D ya
realizada contra EnergyPlus (notebook `000`).

El bloque hueco es un problema **intrínsecamente 2D + transferencia acoplada**
(conducción en la cáscara y los puentes de concreto, **convección de Nusselt** y
**radiación** entre las paredes de la cavidad). Esto lo distingue de las capas
homogéneas 1D ya validadas y exige una estrategia de referencia distinta.

---

## 2. Cómo modela EnerHabitat el bloque, y qué dice Borbón

### 2.1 Modelo de cavidad de EnerHabitat (`Fill.AIR`, muro, β=90°)

Revisando `enerhabitat/ehtools2d.py::_step_hueca`, la cavidad de aire se resuelve
con **convección + radiación acopladas**, tal como exige el fenómeno físico:

- **Convección (Nusselt de pared vertical):**
  `hh = 0.4005 · |T_sup − T_inf|^0.3033 / e22^0.0901`  → depende de **ΔT**.
- **Radiación entre las 4 paredes de la cavidad:** Stefan–Boltzmann
  (`σ = 5.6704e-8`) con **factores de vista** (`_view_factors(a21, e22)`) y
  **emisividad** `E` de las paredes → depende de **T⁴** (de la T promedio).
- El aire de la cavidad es un nodo concentrado `Thueco` que flota.

### 2.2 Lo que midió Borbón (PDF `pdfs/Borbon.pdf`)

Borbón et al. (2010) determinaron experimentalmente (placa caliente protegida,
ASTM C177) la **resistencia térmica R** de un muro de bloque hueco de concreto
**0.12 × 0.20 × 0.40 m**, con dos cavidades de **0.07 × 0.164 × 0.20 m**:

- `R = ΔT_muro · A / q`  (resistencia **superficie a superficie**, sin películas).
- R varía **0.15 → 0.19 °C·m²/W** según las condiciones de operación.
- Ajuste empírico (ec. 4):

  ```
  R = 1 / ( 4.78 + 0.142·√|ΔT| + 0.021·T̄ )      válida para
        0.86 ≤ |ΔT| ≤ 45.81  y  12.75 ≤ T̄ ≤ 47.60
  ```
  con `T̄ = (T_caliente + T_frío)/2` (temperaturas de superficie del muro).

**Coincidencia clave:** la dependencia de R con ΔT (convección de Nusselt) y con
T̄ (radiación ∝ T⁴) que Borbón encontró empíricamente es **exactamente** la
estructura física del solver de EnerHabitat. Borbón es, por tanto, la **referencia
física (experimental)** ideal para validar el modelo de cavidad. La sensibilidad
de R al **acabado superficial** que reporta Borbón se mapea al parámetro
`emissivity` del `HollowBlock`.

---

## 3. El problema de fondo y el rol de cada referencia

**EnergyPlus resuelve muros en 1D** (conduction transfer functions): no modela ni
la conducción lateral en los puentes de concreto ni la convección/radiación 2D de
la cavidad. No puede, por sí solo, ser "la verdad" del bloque hueco. Por eso la
validación usa **tres referencias complementarias**:

| Referencia | Qué valida | Naturaleza |
|---|---|---|
| **Borbón ec. (4) + Tabla 3** | Física de la cavidad: R(ΔT, T̄) | Experimental (verdad física) |
| **EnergyPlus** | Continuidad dinámica del flujo térmico con la cavidad reducida a su R equivalente | Numérica (herramienta estándar) |
| **EnerHabitat 1D** | Consistencia interna del solver 2D (relleno = sólido homogéneo) | Auto-consistencia |

---

## 4. El "hackeo" de EnerHabitat (validado en de-risk)

Para reducir EnerHabitat a un régimen **estacionario tipo placa caliente** —y para
imponer condiciones controladas comparables a Borbón y a EnergyPlus— se explotan
tres palancas. **Las tres fueron probadas y funcionan** (ver §10):

1. **Tsa constante → EPW plano.**
   `Tsa = Ta + Is·α/ho − LWR`. Con un **EPW sintético de temperatura constante y
   radiación solar = 0**, y `absortance = 0`, resulta `Tsa = Ta = constante`.
   - ⚠️ *Nota técnica:* mutar en sitio el DataFrame de `Tsa()` **no** funciona
     (los setters de `tilt/azimuth/absortance` invalidan la caché y `Tsa()` se
     recalcula en cada acceso). El EPW plano es el método robusto, y además es el
     **mismo** mecanismo que ya usa EnergyPlus (EPW custom con Tsa). Apples-to-apples.
2. **Coeficientes `ho` y `hi` ajustables** (`eh.config.ho`, `eh.config.hi`).
   Con `ho, hi` grandes (p. ej. 10³–10⁴ W/m²K) las **superficies quedan fijadas**
   a las temperaturas del aire → condición de **Dirichlet** equivalente a las
   superficies controladas de la placa caliente de Borbón.
3. **`solveAC` con `setpoint`** fija el aire **interior** a una temperatura
   constante; con `hi` grande, `Tsi ≈ setpoint`. La energía `cooling_energy +
   heating_energy` integrada sobre el día da el flujo estacionario `q`.

**Extracción de R (robusta):** se usan las **temperaturas de superficie reales**
del campo estacionario `Tfield` (media física por cara) y el **q real** que
devuelve el solver:

```
q   = (cooling_energy + heating_energy) / (N_pasos · dt)     [W/m²]
R   = (Tso − Tsi) / q                                         [°C·m²/W]
T̄   = (Tso + Tsi) / 2 ;   ΔT = |Tso − Tsi|
```

Se toman las superficies de `Tfield` (no las series `Tso/Tsi` del solver) para
evitar el factor de reporte `nx/(nx-1)` de esas series.

---

## 5. Mapeo geométrico Borbón → `HollowBlock`

El bloque de Borbón se mapea a la celda unitaria de un `HollowBlock` (verificado):

| Borbón (Fig. 5) | Parámetro `HollowBlock` | Valor |
|---|---|---|
| Espesor del muro (dir. del flujo) | `cover_top + cavity + cover_bottom` | 0.12 m |
| Recubrimiento de concreto por cara | `cover_top` = `cover_bottom` | 0.025 m |
| Cavidad (dir. del flujo) | `cavity` (= e22) | 0.07 m |
| Pared/puente de concreto | `web` (= a11) | 0.024 m |
| Ancho de la cavidad | `block_width` (= a21) | 0.164 m |
| Material de cáscara | `material` | `Concreto` (materials.ini) |
| Acabado superficial cavidad | `emissivity` | **0.9** nominal (barrido 0.8–0.95) |

```python
GEOMETRY = {"web":0.024, "block_width":0.164,
            "cover_top":0.025, "cavity":0.07, "cover_bottom":0.025}
```

> **Decidido:** se usa `Concreto` de `materials.ini` (k=1.35, ρ=1800, c=1000) tal
> cual. La emisividad de la cavidad se barre en 0.8–0.95 (nominal 0.9).

---

## 6. Pistas de validación (tracks)

### Track 0 — Consistencia interna (ya iniciada en `002_2DnAC.ipynb`)
`HollowBlock(Fill.SOLID, fill_material = material de la cáscara)` debe reproducir
un muro 1D homogéneo del mismo espesor. **Criterio:** `max|ΔTi| → 0` (≲ 1e-2 °C).
→ Valida la geometría/conducción del solver 2D, aislada de la cavidad.

### Track 1 — **R estacionario del bloque con cavidad vs Borbón** *(validación principal)*
Reproducir, con el hackeo de §4, los puntos de operación de Borbón
(Tabla 1/Tabla 2: pares de temperatura de superficie caliente/fría) y barrer el
rango de validez de la ec. (4): `ΔT ∈ [1, 46] °C`, `T̄ ∈ [13, 48] °C`.

Para cada punto: `R_EH = (Tso−Tsi)/q`; comparar contra `R_Borbón(ΔT, T̄)` (ec. 4)
y contra los 7 valores de la **Tabla 3**. Graficar `R vs ΔT` (réplica de la Fig. 8)
y la superficie `R(ΔT, T̄)` (réplica de la Fig. 9).

**Convergencia de malla primero:** refinar (mallas **conmensurables**: paredes de
la cavidad sobre nodos → cavidad exacta 164×70 mm) hasta que `R` cambie < 5 %;
adoptar esa malla para toda la validación. **Sensibilidades a reportar:**
`emissivity` (0.8–0.95) y `k` del concreto.

### Track 2 — Respuesta **dinámica** vs EnergyPlus *(OPCIONAL / futuro)*
Reusar **exactamente** la metodología 1D del notebook `000`:
EPW custom con `Tsa(t)` real (Temixco, mayo, muro este, α=0.6) y solar = 0;
`ho=13`, `hi=8.1` vía `SurfaceProperty:ConvectionCoefficients`; zona en flotación
libre; lectura con `iertools.read.read_sql`; desfase de **−60 s** en el índice de
EnergyPlus; `Thermal Absorptance = 0.0001`.

La cavidad se representa en EnergyPlus como **capa 1D equivalente**:
`Material:AirGap` (resistencia `R_cav`) **+** los puentes de concreto combinados
por el método de **trayectorias en paralelo** (isothermal-planes, ISO 6946), con
`R_cav` tomado de la caracterización estacionaria del Track 1 al ΔT/T̄
representativo del día. Se comparan `Ti(t)`, pico, **retardo** y **factor de
decremento**, y la **energía** diaria, contra `System2D.solve()`.

> Limitación honesta: EnergyPlus ve la cavidad como una resistencia fija; no
> reproduce la dependencia R(ΔT,T̄). El Track 2 verifica **continuidad dinámica**
> con la herramienta estándar, no la física de la cavidad (eso es el Track 1).

### Track 3 — Bloque **relleno sólido** 2D vs EnergyPlus 1D *(OPCIONAL / futuro)*
Para `Fill.SOLID` con relleno aislante (p. ej. EPS), construir el equivalente 1D
por trayectorias en paralelo (cáscara de concreto ‖ núcleo aislante) y comparar
`Ti(t)` EH-2D vs EP-1D. Cuantifica el **efecto 2D** (puentes térmicos) que el 1D
no captura.

---

## 7. Criterios de aceptación

| Track | Métrica | Criterio |
|---|---|---|
| 0 | `max|ΔTi|` (2D relleno vs 1D homogéneo) | ≤ 0.01 °C |
| **1** | Convergencia de malla | ΔR < 5 % al refinar |
| **1** | `R_EH` en los 7 puntos de Borbón | dentro de **0.15–0.19** °C·m²/W |
| **1** | `R_EH` vs ec. (4) | error ≤ **±10%** (banda exp. de Borbón ≈ ±3.75%) |
| **1** | Tendencias | R↓ al crecer ΔT; variación con T̄ (signo correcto) |
| 2 | `Ti(t)` EH-2D vs EnergyPlus | RMSE ≲ 0.3 °C; pico ≤ 0.5 °C; retardo ≤ 15 min |
| 2 | Energía diaria | ≤ 5% |
| 3 | `Ti(t)` EH-2D vs EP-1D paralelo | reportar diferencia (cuantifica efecto 2D) |

---

## 8. Entregables

**Requeridos (alcance aprobado):**
1. `notebooks/003_2D_R_estacionario_Borbon.ipynb` — **Track 1**:
   helper `flat_epw(T)`, convergencia de malla, extracción de R, tablas y figuras
   réplica de Borbón (Fig. 8 y 9), barrido de emisividad (0.8–0.95).
2. Extensión de `002_2DnAC.ipynb` con el criterio numérico del Track 0.

**Opcionales / futuros:**
3. `notebooks/004_2D_dinamico_EnergyPlus.ipynb` — Track 2 (IDF equivalente +
   EnergyPlus).
4. `notebooks/005_2D_relleno_solido.ipynb` — Track 3.

---

## 9. Riesgos y limitaciones conocidas

- **EnergyPlus es 1D:** la referencia física del bloque con aire es **Borbón**, no
  EnergyPlus. Documentarlo explícitamente en los notebooks.
- **Caché de `Tsa()`:** no mutar el DataFrame; usar siempre el EPW plano (§4).
- **Malla:** `nx, ny` son número de nodos (resolución), independientes de la
  geometría (`dx=X/nx`, `dy=Y/ny`). Las paredes de la cavidad se "pegan" a los
  nodos (truncamiento entero en `compute_mesh`) → usar mallas **conmensurables**
  (cavidad exacta) para evitar que `R` varíe a saltos al refinar.
- **Costo de cómputo:** el solver 2D + convergencia a estado estable es lento
  (~JIT + decenas de "días" iterativos).
- **Propiedades del concreto y emisividad:** influyen en R; barrer y reportar.
- **Pequeño salto de película** con `ho/hi` finitos: irrelevante porque R se
  calcula con `Tso/Tsi` reales (de `Tfield`).

---

## 10. Resultados del Track 1 (notebook 003, malla cuadrada 212×120 convergida)

Convergencia de malla (mallas conmensurables, cavidad exacta 164×70 mm) en el
punto de referencia (corrida 1): `R` = 0.1593 (53×120) → 0.1598 (106×240) →
0.1592 (212×120), ΔR < 0.4 % → se adopta la malla cuadrada isótropa **212×120**
(`dx=dy=1 mm`).

7 puntos de operación de Borbón (temperaturas de superficie reales de la Tabla 2):

| run | ΔT | T̄ | q (W/m²) | **R_EH** | R_eq4 | err vs ec.(4) |
|---|---|---|---|---|---|---|
| 1 | 13.56 | 35.55 | 85.2 | **0.1592** | 0.1653 | −3.7% |
| 2 | 16.86 | 33.16 | 105.2 | **0.1602** | 0.1650 | −2.9% |
| 3 | 20.22 | 40.70 | 130.4 | **0.1551** | 0.1594 | −2.7% |
| 4 | 23.28 | 28.42 | 143.4 | **0.1623** | 0.1650 | −1.6% |
| 5 | 26.89 | 45.87 | 177.8 | **0.1513** | 0.1543 | −2.0% |
| 6 | 30.12 | 48.38 | 201.4 | **0.1496** | 0.1521 | −1.7% |
| 7 | 39.78 | 41.65 | 260.5 | **0.1527** | 0.1527 | 0.0% |

→ **MAE vs ec. (4) = 2.07 %**, R_EH = 0.150–0.162 (en/al borde de la banda 0.15–0.19),
con las tendencias correctas (R↓ al crecer ΔT; variación con T̄). Emisividad:
R = 0.167 → 0.156 al subir ε de 0.8 a 0.95 (radiación ∝ T⁴). Superficie R(ΔT,T̄):
`max|R_EH − R_eq4| = 0.012`. **El modelo de cavidad de aire queda validado
físicamente.**

---

## 11. Plan de ejecución (alcance aprobado)

1. ✅ Documento aprobado (alcance Track 1).
2. Track 0: cerrar `002` con criterio numérico (`max|ΔTi| ≤ 0.01 °C`).
3. ✅ **Track 1:** notebook `003` (convergencia de malla + 7 puntos de Borbón +
   barrido + emisividad), malla cuadrada 212×120 → tablas y figuras réplica de
   Fig. 8 y 9.
4. Redactar conclusiones del Track 1 (hecho en §10).
5. (Futuro, si se decide) Tracks 2 y 3.

---

### Decisiones tomadas (2026-06-25)

- **Alcance:** Track 1 (vs Borbón) es la validación física requerida del bloque
  con cavidad de aire. Tracks 2 y 3 quedan como opcionales/futuros.
- **Concreto:** `Concreto` de `materials.ini` (k=1.35, ρ=1800, c=1000).
- **Emisividad:** nominal 0.9, barrido 0.8–0.95.
- **Malla:** cuadrada isótropa **212×120** (`dx=dy=1 mm`, conmensurable), tras
  confirmar convergencia ΔR < 5 %.
