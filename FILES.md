# Tabla de control de archivos — Célula Virtual

Documento vivo. Se actualiza cada vez que se crea o modifica un archivo relevante del proyecto. Referencia de diseño: `whitepaper_celula_virtual.md`.

**Estado**: 🔲 Pendiente · 🟡 En progreso · ✅ Implementado y validado

| Archivo | Funcionalidad principal | Capa | Fase | Estado |
|---|---|---|---|---|
| `engine/species.py` | Gestor de especies moleculares (catálogo y cantidades) | 1 — Motor | 0 | ✅ |
| `engine/reactions.py` | Gestor de reacciones (reactivos, productos, propensidad) | 1 — Motor | 0 | ✅ |
| `engine/ssa.py` | Bucle SSA propio (algoritmo de Gillespie, generador/iterador) | 1 — Motor | 0 | ✅ |
| `engine/trajectory.py` | Registro de trayectorias temporales (pandas) | 1 — Motor | 0 | ✅ |
| `engine/debugger.py` | Ejecución paso a paso, breakpoints condicionales, inspección de estado | 1 — Motor | 0.5 | ✅ |
| `engine/propensity_view.py` | Visualización en tiempo real del espacio de propensidades | 1 — Motor | 5.1 | 🔲 |
| `engine/forensics.py` | Análisis causal retrospectivo ("modo forense", sección 4.5) | 1 — Motor | 5.2 | 🔲 |
| `engine/snapshot.py` | Snapshots de estado + retroceso (undo) | 1 — Motor | 5.3 | 🔲 |
| `biology/genome.py` | `GenomeModule` — instala especies/reacciones de transcripción (simplificada, un paso) a partir de un genoma cargado | 2 — Biología | 1 | ✅ |
| `biology/transcription.py` | `TranscriptionModule` — ARN-pol, elongación, terminación | 2 — Biología | 1 | 🔲 |
| `biology/translation.py` | `TranslationModule` — instala especies de proteína y reacciones de traducción (simplificada, un paso) sobre genes ya transcritos | 2 — Biología | 1 | ✅ |
| `biology/degradation.py` | `DegradationModule` — degradación de ARNm y proteínas | 2 — Biología | 1 | 🔲 |
| `biology/metabolism.py` | `MetabolismModule` — glucólisis simplificada, acoplamiento energético | 2 — Biología | 3 | 🔲 |
| `biology/replication.py` | `ReplicationModule` — replicación de ADN y división celular | 2 — Biología | 4 | 🔲 |
| `data_io/genome_loader.py` | Carga de genoma desde FASTA + anotación GFF3/CSV (sección 6.2) | I/O | 0 | ✅ |
| `data_io/sbml_io.py` | Import/export de modelos en formato SBML (sección 7) | I/O | — | 🔲 |
| `cell.py` | Orquestador `Cell` — compone módulos, API pública, modo normal y depurado | 3 — Orquestador | 0 / 0.5 | ✅ (`run()`, `get_trajectory()`, `debug()`) |
| `tests/test_ssa.py` | Validación mecánica (conservación, monotonía, reproducibilidad) y estadística (tiempo de espera exponencial, elección proporcional a propensidad) del motor SSA — 10 casos | Tests | 0 | ✅ |
| `tests/test_species.py` | Validación unitaria de `engine/species.py` (13 casos, sin pytest disponible en este entorno: ejecutados con `assert` puro) | Tests | 0 | ✅ |
| `tests/test_reactions.py` | Validación unitaria de `engine/reactions.py` — propensidades uni/bimoleculares, `apply()`, validaciones (13 casos) | Tests | 0 | ✅ |
| `tests/test_cell.py` | Validación unitaria de `cell.py` — encadenado, `run()`, `max_steps`/`max_time`, `get_trajectory()`, `debug()` y su sincronización con `run()` (13 casos) | Tests | 0 / 0.5 | ✅ |
| `tests/test_trajectory.py` | Validación de `engine/trajectory.py` — reconstrucción, conservación, columnas, errores (7 casos) | Tests | 0 | ✅ |
| `tests/test_debugger.py` | Validación de `engine/debugger.py` — incluye el criterio central: paso a paso produce exactamente la misma secuencia que `run()` con el mismo seed (13 casos) | Tests | 0.5 | ✅ |
| `tests/test_genome_loader.py` | Validación de `data_io/genome_loader.py` — carga, hebra +/-, y todas las validaciones de FASTA/anotación (14 casos) | Tests | 0 | ✅ |
| `tests/test_genome.py` | Validación de `biology/genome.py` — instalación, conflictos de nombre, tasas por promotor, integración de extremo a extremo con el genoma de juguete real (10 casos) | Tests | 1 | ✅ |
| `tests/test_translation.py` | Validación de `biology/translation.py` — instalación, conflictos, tasas por RBS, pipeline completo ADN→ARNm→proteína (10 casos) | Tests | 1 | ✅ |
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
