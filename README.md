# Célula Virtual: Simulación Biofísica de una Célula Mínima con Depuración Interactiva

Simulación biofísica de una célula mínima, con motor de cinética estocástica propio y un **depurador interactivo** — pausa, breakpoints, inspección del espacio de propensidades, análisis causal retrospectivo y retroceso (undo) — como diferenciador frente a otros whole-cell models existentes.

> Proyecto en desarrollo activo. El diseño completo, las decisiones tomadas y el porqué de cada una viven en [`docs/whitepaper_celula_virtual.md`](docs/whitepaper_celula_virtual.md) — es el documento de referencia de este repositorio, no solo este README.

## Qué es esto

Un software que representa una célula virtual vacía, a la que se le puede añadir una secuencia de ADN real (FASTA + anotación) y observar su evolución temporal — expresión génica, regulación, metabolismo acoplado a energía, y división celular — simulada con el algoritmo de Gillespie (cinética química estocástica exacta), no con reglas inventadas.

Lo que lo diferencia de otros simuladores de biología de sistemas (E-Cell, VCell, los whole-cell models de Karr/Covert Lab) es que el motor de simulación se escribió desde cero, específicamente para poder **inspeccionarlo por dentro mientras corre**, en vez de tratarlo como una caja negra que solo se ejecuta y produce una gráfica al final.

## Funcionalidades

- **Motor de simulación propio** (`engine/`): algoritmo de Gillespie exacto, escrito como iterador para poder pausarse en cualquier punto.
- **Depurador interactivo** (`engine/debugger.py`): ejecución paso a paso, breakpoints condicionales sobre el estado molecular, visualización en vivo del espacio de propensidades, análisis causal retrospectivo ("¿qué consumió todo el ATP?"), y retroceso (undo) reproducible bit a bit.
- **Expresión génica** (`biology/`): transcripción, traducción, degradación, regulación (represión), metabolismo acoplado a ATP, y división celular con reparto binomial de especies entre células hijas.
- **Carga de genomas reales** (`data_io/genome_loader.py`): FASTA + anotación CSV — nunca un genoma codificado a mano — con diseño pensado para escalar de 5 a cientos de genes sin cambiar código.
- **Interoperabilidad SBML** (`data_io/sbml_io.py`): exportar/importar modelos en formato estándar.
- **Interfaz gráfica de escritorio** (`gui/`, PySide6, tema oscuro): todas las funcionalidades anteriores accesibles desde una ventana con pestañas — Genoma, Ejecución, Depurador, Forense.

## Instalación

```bash
git clone <url-de-este-repositorio>
cd celula-virtual
pip install -r requirements.txt
```

Requiere Python 3.10+. Ver [`requirements.txt`](requirements.txt) para el detalle de qué librería hace falta en cada fase del proyecto.

## Uso rápido

### Interfaz gráfica

```bash
python -m gui.app
```

Desde la pestaña **Genoma**, carga `examples/genomas/toy_genome.fasta` y su anotación, instala, y pasa a **Ejecución** para correr la simulación.

### Desde código

```python
from cell import Cell
from engine.reactions import Reaction

cell = Cell(seed=1)
cell.add_species("A", 100).add_species("B", 0)
cell.add_reaction(Reaction("A_to_B", {"A": 1}, {"B": 1}, rate_constant=1.0))

cell.run()
print(cell.species.get_count("B"))  # 100

df = cell.get_trajectory()  # pandas.DataFrame, tiempo x especie x cantidad
```

### Depurador

```python
debugger = cell.debug(enable_snapshots=True)
debugger.set_breakpoint("poca_A", lambda s: s["A"] < 10)
debugger.run_until_breakpoint()
print(debugger.pending_propensities())  # qué reacciones podían haber disparado
debugger.undo_to(5)  # retroceder al paso 5, exactamente
```

### Un genoma real

```python
from data_io.genome_loader import load_genome
from biology.genome import GenomeModule
from biology.translation import TranslationModule
from biology.degradation import DegradationModule

genome = load_genome("examples/genomas/toy_genome.fasta", "examples/genomas/toy_genome_annotation.csv")
cell = Cell(seed=1)
installed = GenomeModule(genome).install(cell)
translated = TranslationModule(installed).install(cell)
DegradationModule(installed, translated).install(cell)  # sin esto, ARNm/proteína crecen sin límite

cell.run(max_steps=500)
```

## Ejecutar los tests

```bash
pytest tests/ -v
```

## Estructura del proyecto

```
celula-virtual/
├── engine/        # Motor de simulación (SSA, depurador, snapshots, forense, trayectorias)
├── biology/       # Módulos biológicos (genoma, transcripción, traducción, degradación, regulación, metabolismo, replicación)
├── data_io/       # Carga de genomas (FASTA + anotación) e interoperabilidad SBML
├── gui/           # Interfaz gráfica de escritorio (PySide6)
├── tools/         # Utilidades de desarrollo (generador de genomas sintéticos)
├── tests/         # Suite de tests
├── examples/      # Genomas de ejemplo
├── docs/          # Whitepaper de diseño
├── requirements.txt
└── FILES.md       # Mapa vivo de qué archivo hace qué y su estado
```

Para saber exactamente qué hace cada archivo y en qué estado está (implementado, en progreso, pendiente), consulta [`FILES.md`](FILES.md) — es un documento vivo que se actualiza en cada cambio relevante.

## Documentación

- [`docs/whitepaper_celula_virtual.md`](docs/whitepaper_celula_virtual.md) — diseño completo: motivación, arquitectura, decisiones y por qué se tomaron, roadmap por fases, limitaciones conocidas.
- [`FILES.md`](FILES.md) — tabla de control de archivos.

## Estado del proyecto

El ciclo celular central (motor SSA, depurador completo, expresión génica, regulación, metabolismo, división celular) y la interoperabilidad SBML están implementados y validados con tests. La interfaz gráfica está construida pero pendiente de validación funcional completa en un entorno con pantalla real (ver `FILES.md`). El roadmap detallado, fase a fase, está en la sección 8 del whitepaper.

## Limitaciones conocidas

Documentadas en detalle en la sección 12 del whitepaper. Las más relevantes:

- Simulación no espacial (sin difusión ni geometría celular real).
- Sin mutación del ADN — dinámica temporal, no evolución darwiniana (decisión deliberada, no limitación técnica).
- Transcripción y traducción modeladas como un único paso estocástico, no como la cadena completa de unión/elongación/terminación.
- La GUI no soporta simular varias células a la vez tras una división.

## Licencia

Por decidir (candidatas: MIT o GPLv3 — esta última es la que usa E-Cell4, lo que facilitaría eventuales colaboraciones cruzadas). Ver sección 10 del whitepaper.

## Contribuir

Proyecto en desarrollo activo, todavía sin proceso formal de contribución externa. Si quieres reportar un problema o sugerir algo, abre un issue.
