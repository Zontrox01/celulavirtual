# Tabla de control de archivos — Célula Virtual

Documento vivo. Se actualiza cada vez que se crea o modifica un archivo relevante del proyecto. Referencia de diseño: `whitepaper_celula_virtual.md`.

**Estado**: 🔲 Pendiente · 🟡 En progreso · ✅ Implementado y validado

| Archivo | Funcionalidad principal | Capa | Fase | Estado |
|---|---|---|---|---|
| `engine/species.py` | Gestor de especies moleculares (catálogo y cantidades) | 1 — Motor | 0 | ✅ |
| `engine/reactions.py` | Gestor de reacciones (reactivos, productos, propensidad) | 1 — Motor | 0 | ✅ |
| `engine/ssa.py` | Bucle SSA propio (algoritmo de Gillespie, generador/iterador); `get_full_state()`/`restore_full_state()` para el retroceso | 1 — Motor | 0 / 5.3 | ✅ |
| `engine/trajectory.py` | Registro de trayectorias temporales (pandas) | 1 — Motor | 0 | ✅ |
| `engine/debugger.py` | Ejecución paso a paso, breakpoints condicionales, inspección de estado, `undo_to()` (retroceso opcional vía snapshots) | 1 — Motor | 0.5 / 5.3 | ✅ |
| `engine/propensity_view.py` | Visualización del espacio de propensidades (datos testeables + render `rich` con fallback en texto) | 1 — Motor | 5.1 | ✅ |
| `engine/forensics.py` | Análisis causal retrospectivo ("modo forense", sección 4.5) — `ForensicsAnalyzer.analyze_species()` | 1 — Motor | 5.2 | ✅ |
| `engine/snapshot.py` | Snapshots de estado + retroceso (undo) — `SnapshotStore`, reproducción determinista desde el snapshot más cercano | 1 — Motor | 5.3 | ✅ |
| `biology/genome.py` | `GenomeModule` — instala especies/reacciones de transcripción; acoplamiento energético opcional a ATP (Fase 3) | 2 — Biología | 1 / 3 | ✅ |
| `biology/transcription.py` | `TranscriptionModule` — ARN-pol, elongación, terminación | 2 — Biología | 1 | 🔲 |
| `biology/translation.py` | `TranslationModule` — instala proteínas y reacciones de traducción; acoplamiento energético opcional a ATP (Fase 3) | 2 — Biología | 1 / 3 | ✅ |
| `biology/degradation.py` | `DegradationModule` — reacciones de degradación de un paso para ARNm y proteínas, por vida media | 2 — Biología | 1 | ✅ |
| `biology/regulation.py` | `RegulationModule` — represión transcripcional reversible (unión/disociación represor-gen), no prevista en el diseño original | 2 — Biología | 2 | ✅ |
| `biology/metabolism.py` | `MetabolismModule` — importación de glucosa + glicólisis simplificada (opcionalmente catalizada por una enzima ya traducida), produce el ATP del acoplamiento energético | 2 — Biología | 3 | ✅ |
| `biology/replication.py` | `ReplicationModule` — condición de división + reparto binomial de especies entre 2 células hijas | 2 — Biología | 4 | ✅ |
| `data_io/genome_loader.py` | Carga de genoma desde FASTA + anotación GFF3/CSV (sección 6.2) | I/O | 0 | ✅ |
| `data_io/sbml_io.py` | Import/export de modelos en formato SBML (sección 7) | I/O | — | 🔲 |
| `cell.py` | Orquestador `Cell` — compone módulos, API pública, modo normal y depurado, incluido `undo` sincronizado | 3 — Orquestador | 0 / 0.5 / 5.3 | ✅ (`run()`, `get_trajectory()`, `debug()` con `on_event`/`on_undo` sincronizados) |
| `tests/test_ssa.py` | Validación mecánica (conservación, monotonía, reproducibilidad) y estadística (tiempo de espera exponencial, elección proporcional a propensidad) del motor SSA — 10 casos | Tests | 0 | ✅ |
| `tests/test_species.py` | Validación unitaria de `engine/species.py` (13 casos, sin pytest disponible en este entorno: ejecutados con `assert` puro) | Tests | 0 | ✅ |
| `tests/test_reactions.py` | Validación unitaria de `engine/reactions.py` — propensidades uni/bimoleculares, `apply()`, validaciones (13 casos) | Tests | 0 | ✅ |
| `tests/test_cell.py` | Validación unitaria de `cell.py` — encadenado, `run()`, `max_steps`/`max_time`, `get_trajectory()`, `debug()` y su sincronización con `run()` (13 casos) | Tests | 0 / 0.5 | ✅ |
| `tests/test_trajectory.py` | Validación de `engine/trajectory.py` — reconstrucción, conservación, columnas, errores (7 casos) | Tests | 0 | ✅ |
| `tests/test_debugger.py` | Validación de `engine/debugger.py` — incluye el criterio central: paso a paso produce exactamente la misma secuencia que `run()` con el mismo seed (13 casos) | Tests | 0.5 | ✅ |
| `tests/test_snapshot.py` | Validación de `engine/snapshot.py` — intervalos, snapshot más cercano, y el criterio central: retroceder + reproducir = exactamente el futuro ya ocurrido (10 casos) | Tests | 5.3 | ✅ |
| `tests/test_debugger_undo.py` | Validación de la integración `undo_to()` en `Debugger` — truncado de histórico, sincronización de `has_ended`, reproducibilidad exacta, regresión (9 casos) | Tests | 5.3 | ✅ |
| `tests/test_cell_debug_undo.py` | Cierra la limitación de sincronización Cell↔undo — recorte de `cell.get_events()`, `get_trajectory()` consistente, continuar con `run()` tras un retroceso (7 casos) | Tests | 5.3 | ✅ |
| `tests/test_genome_loader.py` | Validación de `data_io/genome_loader.py` — carga, hebra +/-, y todas las validaciones de FASTA/anotación (14 casos) | Tests | 0 | ✅ |
| `tests/test_genome.py` | Validación de `biology/genome.py` — instalación, conflictos de nombre, tasas por promotor, acoplamiento ATP opcional, integración de extremo a extremo (13 casos) | Tests | 1 / 3 | ✅ |
| `tests/test_translation.py` | Validación de `biology/translation.py` — instalación, conflictos, tasas por RBS, acoplamiento ATP opcional, pipeline ADN→ARNm→proteína (13 casos) | Tests | 1 / 3 | ✅ |
| `tests/test_degradation.py` | Validación de `biology/degradation.py` — tasas por vida media, pipeline de 3 módulos alcanzando equilibrio acotado (10 casos) | Tests | 1 | ✅ |
| `tests/test_regulation.py` | Validación de `biology/regulation.py` — toggle libre/reprimido, y comparación con/sin represión activa (criterio de la Fase 2) (10 casos) | Tests | 2 | ✅ |
| `tests/test_metabolism.py` | Validación de `biology/metabolism.py` — glicólisis catalítica, criterio de la Fase 3 (524 vs 5 ARNm con/sin ATP), breakpoint de ejemplo del whitepaper, ciclo cerrado gen→enzima→metabolismo (14 casos) | Tests | 3 | ✅ |
| `tests/test_replication.py` | Validación de `biology/replication.py` — condición de división, reparto binomial, ciclo celular completo hasta dos hijas independientes (13 casos) | Tests | 4 | ✅ |
| `tests/test_propensity_view.py` | Validación de `engine/propensity_view.py` — preparación de datos, render en texto, fallback real sin `rich` instalado, integración con el depurador (11 casos) | Tests | 5.1 | ✅ |
| `tests/test_forensics.py` | Validación de `engine/forensics.py` — ventanas temporales, contribuciones netas, catalizadores excluidos, identificación real de productor/consumidor en una simulación (12 casos) | Tests | 5.2 | ✅ |
| `examples/genomas/toy_genome.fasta` | Genoma de juguete: 4 genes (2 ribosomales, 1 metabólico, 1 regulador) | Datos | 1 | ✅ |
| `examples/genomas/toy_genome_annotation.csv` | Anotación de los 4 genes (posición, hebra, promoter_id) | Datos | 1 | ✅ |
| `docs/whitepaper_celula_virtual.md` | Documento de diseño del proyecto (decisiones, arquitectura, roadmap) | Docs | — | ✅ |
| `docs/FILES.md` | Este documento — mapa vivo del repositorio | Docs | — | 🟡 |
| `requirements.txt` | Dependencias Python del proyecto, agrupadas por fase (whitepaper, sección 7.1) | Docs | 0 | ✅ |

---

## Cómo mantener esta tabla

- Al crear un archivo nuevo con una responsabilidad ya prevista en el whitepaper (sección 5.4): actualizar su fila de 🔲 a 🟡 (en progreso) y luego a ✅ (implementado y validado, es decir, con sus tests correspondientes pasando).
- Al crear un archivo con una responsabilidad **no prevista** en el diseño original: añadir una fila nueva y valorar si el whitepaper necesita una actualización correspondiente (evitar que el código se desvíe del documento de diseño sin que quede constancia).
- Si un archivo empieza a mezclar más de una responsabilidad (regla de la sección 5.4), se divide y esta tabla se actualiza con las filas resultantes.

---

## Nota: bug corregido en `engine/ssa.py` (detectado durante `biology/degradation.py`)

Al construir el test de equilibrio dinámico de `degradation.py`, se detectó que alcanzar `max_time` marcaba la simulación como terminada de forma **permanente** (`_ended = True`), impidiendo reanudarla aunque se subiera `max_time` después — rompiendo la promesa de `cell.run()` de poder llamarse varias veces seguidas para seguir avanzando la misma simulación.

**Corrección**: se introdujo `PropensityExhausted` (subclase de `SimulationEnded`) para el caso realmente permanente (propensidad total = 0), dejando `SimulationEnded` genérico para el caso de `max_time`, que ahora es una pausa reanudable. `engine/debugger.py` se actualizó igual, distinguiendo ambos casos en `has_ended`.

Tests de regresión añadidos: `test_max_time_is_a_resumable_pause_not_a_permanent_end` y `test_exhaustion_remains_permanent_unlike_max_time` en `test_ssa.py`; `test_max_time_pause_does_not_mark_debugger_as_permanently_ended` en `test_debugger.py`.

---

## Nota: bug corregido en `cell.py` (detectado al cerrar la sincronización Cell↔undo)

Al implementar `_truncate_events_to()` para que `cell.debug().undo_to()` recortara también el histórico de la propia `Cell`, la primera versión **reasignaba** `self._events` a una lista nueva en cada recorte. Pero el callback `on_event` que `Cell.debug()` le pasa al `Debugger` es `self._events.append` — un método ya vinculado a la lista *original* en el momento de crear el `Debugger`. Al reasignar `self._events` a una lista nueva, ese callback seguía apuntando a la lista vieja, huérfana: cualquier evento nuevo tras un `undo_to()` dejaba de aparecer en `cell.get_events()`, aunque sí aparecía en `debugger.history()`.

**Corrección**: `_truncate_events_to()` muta la lista en el mismo sitio (`self._events[:] = [...]`, slice assignment) en vez de reasignarla, preservando la referencia que el callback ya tiene vinculada.

Test de regresión: `test_multiple_undo_and_redo_cycles_stay_consistent` en `test_cell_debug_undo.py` — sin este test (que encadena varios ciclos de avance/retroceso) el bug no se habría detectado, porque un solo retroceso sin avanzar después no lo manifiesta.
