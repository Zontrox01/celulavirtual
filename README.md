# celulavirtual
Célula Virtual: Simulación Biofísica de una Célula Mínima con Depuración Interactiva

White paper de diseño — v0.2 Fecha: Agosto 2026

Resumen ejecutivo

Este documento define el diseño de un software en Python capaz de representar una célula virtual vacía, a la que se le puede añadir una secuencia de ADN, y observar la dinámica temporal resultante: expresión génica, metabolismo energético, replicación del ADN y división celular, simulado con cinética química estocástica (algoritmo de Gillespie).

El proyecto se enmarca en la tradición de los whole-cell models de biología de sistemas — en particular el modelo de Mycoplasma genitalium (Karr et al., 2012) y su sucesor sobre la célula mínima sintética JCVI-syn3A (Covert Lab, Stanford, 2022) — pero incorpora un elemento de innovación no presente en las herramientas existentes: un depurador interactivo para la simulación biológica, que permite pausar la ejecución, inspeccionar reacción a reacción, y establecer condiciones de parada ("breakpoints") sobre el estado molecular de la célula — trayendo la metáfora de la depuración de software al terreno de la biología de sistemas.

No es un simulador de evolución darwiniana. No hay mutación ni selección de poblaciones. El objetivo es observar y poder inspeccionar en detalle el ciclo de vida de una célula según las leyes de la cinética química.

1. Motivación y precedentes
1.1 El estado del arte

Los modelos de célula completa más relevantes son:

E-Cell (desde 1996) y su versión actual E-Cell4: plataforma de código abierto para modelar sistemas multiescala como la célula, con núcleo en C++ y frontend en Python.
Virtual Cell (VCell): software libre con interfaz gráfica sin necesidad de programar; el usuario introduce reacciones y VCell genera el modelo matemático automáticamente. Permite ejecutar en local o delegar a sus servidores.
Whole-Cell Model de M. genitalium (Karr et al., Cell, 2012): primer modelo que integró los procesos biológicos conocidos de un organismo completo (525 genes).
Modelo de JCVI-syn3A (Covert Lab, Stanford, 2022–2023): sucesor sobre la célula mínima sintética (~473 genes), extendido posteriormente a colonias de E. coli mediante la librería Vivarium, que permite simulaciones paralelizadas en la nube con detalle molecular completo en un entorno espacial compartido.
1.2 El hueco identificado

Revisando la literatura sobre los obstáculos actuales del campo, aparecen dos problemas recurrentes:

La construcción fiable de los modelos y la inferencia de parámetros son señaladas como uno de los mayores cuellos de botella del campo; se necesitan métodos más eficientes para simular y depurar sistemas bioquímicos estocásticos grandes.
La propia comunidad de whole-cell modeling pide explícitamente mejores herramientas de control de versiones, validación y trazabilidad para estos modelos, algo que hoy no existe de forma satisfactoria.

Ninguna de las herramientas existentes (E-Cell, VCell, los modelos de Karr/Covert Lab) ofrece una forma de inspeccionar el mecanismo interno de la simulación paso a paso, de manera análoga a como un desarrollador depura un programa. Todas presentan la simulación como una caja relativamente cerrada: se define el modelo, se ejecuta, y se analiza el resultado final en forma de gráficas de trayectorias. Ahí está el hueco que este proyecto quiere ocupar.

2. Objetivos y alcance
2.1 Objetivo general

Construir un motor de simulación en Python que modele el comportamiento temporal de una célula mínima a partir de una secuencia de ADN proporcionada por el usuario, con base biofísica rigurosa, y que incorpore un depurador interactivo como forma primaria de explorar e interpretar la simulación.

2.2 Objetivos específicos (MVP)
Motor de simulación estocástica (Gillespie SSA) escrito desde cero, diseñado explícitamente para exponer su estado interno paso a paso.
Depurador interactivo: ejecución paso a paso, breakpoints sobre condiciones del estado molecular, inspección de propensidades de todas las reacciones candidatas en cada instante.
Módulo de expresión génica: transcripción y traducción realistas a partir de una secuencia ADN real (codones, marco de lectura, código genético estándar).
Módulo metabólico mínimo que provea energía (ATP) y precursores (NTPs, aminoácidos) a los procesos anteriores.
Módulo de replicación de ADN y división celular con segregación estocástica de moléculas.
Registro y visualización de trayectorias temporales, enlazado con el histórico de eventos del depurador.
Diseño escalable del genoma de entrada: el número de genes debe ser un parámetro de datos (archivo de secuencia + anotación), no una constante del código, para que el mismo software sirva después con genomas de cientos de genes en estudios más serios.
2.3 Explícitamente fuera de alcance (por ahora)
Simulación espacial (difusión 3D, geometría celular real).
Mutación del ADN y evolución darwiniana / poblacional.
Redes metabólicas a escala genómica completa (FBA de miles de reacciones).
Interacciones célula-célula o entornos multicelulares.
Construcción de modelos asistida por LLM (considerada, aparcada para una fase futura — ver sección 9).
Control de versiones tipo "git" para modelos biológicos (considerada, aparcada para una fase futura — ver sección 9).
3. Decisiones de diseño fundamentales
Decisión	Elección	Justificación
Espacialidad	No espacial (well-mixed)	Reduce drásticamente la complejidad; primera aproximación válida y extensible después.
Naturaleza de la cinética	Estocástica (Gillespie SSA)	Los números de copias reales de ARNm y reguladores son bajos; las ODEs deterministas no capturan el ruido biológico, que es parte esencial del fenómeno.
Alcance temporal	Dinámica de una célula (o linaje), sin mutación	El interés es observar el ciclo de vida, no la evolución darwiniana.
Implementación del motor SSA	Propia, desde cero, en Python	Librerías maduras como GillesPy2 ejecutan la simulación como una caja negra optimizada en C++/Cython, invocada mediante un único método .run(). Esto las hace rápidas pero imposibles de pausar o inspeccionar a mitad de ejecución sin modificar su código interno. Como el depurador es el diferenciador del proyecto, el bucle SSA debe ser nuestro desde el principio: es donde vive el punto de enganche.
Resto de componentes no críticos para la innovación	Librerías de terceros maduras	Reinventar el manejo de secuencias de ADN, formatos de intercambio o graficado no aporta valor innovador — solo introduce riesgo de errores evitables.
Fuente del "ADN"	Secuencia real (BioPython Seq), con anotación de genes	Evita inventar reglas arbitrarias; usa el código genético estándar real.
Entrada del genoma	Archivo externo (FASTA + anotación GFF3/CSV), nunca codificado en Python	Permite escalar de 8 a cientos de genes sin tocar código, y cargar genomas reales descargados de bases de datos públicas.
Formato de red de reacciones	Compatible con SBML donde sea posible (import/export)	Permite interoperar con herramientas existentes y comparar contra modelos publicados.
4. El diferenciador: depurador interactivo
4.1 Concepto

Igual que un IDE permite pausar un programa, poner un breakpoint, e inspeccionar variables en memoria, este depurador permite:

Ejecución paso a paso: avanzar reacción a reacción (no solo intervalo de tiempo a intervalo de tiempo), viendo exactamente qué reacción disparó, con qué probabilidad relativa frente a las demás candidatas, y qué especies cambiaron.
Breakpoints condicionales: pausar automáticamente cuando se cumple una condición sobre el estado molecular — por ejemplo, ATP < 50, o proteína_X == 0 y ARNm_X > 0 (indicio de que la traducción se ha detenido pese a haber ARN mensajero disponible).
Inspección del estado completo: en cualquier pausa, ver el vector completo de especies y sus cantidades, así como las propensidades de todas las reacciones que podrían haber ocurrido en ese instante (no solo la que ocurrió).
Historial navegable: retroceder sobre el log de eventos ya ocurridos, no solo avanzar.
4.2 Por qué es viable con recursos modestos

El coste computacional de esta funcionalidad es marginal: la información que expone el depurador (qué reacción disparó, con qué propensidad, qué cambió) ya se calcula internamente en cada iteración de cualquier implementación correcta del algoritmo de Gillespie. El depurador no añade cálculo nuevo — añade una capa de control (pausar/reanudar) y de registro (loggear el estado en cada paso) sobre un bucle que de todos modos hay que ejecutar. En un genoma de juguete de 5-8 genes, con pocas docenas de reacciones activas, esto es trivial para un equipo con un Ryzen 5 Pro y 8 GB de RAM.

4.3 Relación con la interpretabilidad futura

El log de eventos que genera el depurador (secuencia ordenada de: instante, reacción disparada, estado resultante) es exactamente la materia prima que necesitaría, en una fase futura, un sistema de "narrativa causal" que traduzca la traza estocástica a explicaciones en lenguaje natural (ver sección 9). Diseñar el depurador pensando en esta reutilización evita tener que rehacer el registro de eventos más adelante.

4.4 Funcionalidades avanzadas del depurador

Tres capacidades adicionales que profundizan el diferenciador del proyecto:

a) Visualización en tiempo real del "espacio de propensidades" En cada paso, el motor SSA calcula la propensidad de todas las reacciones candidatas, no solo de la que finalmente dispara. Hoy esa información se descarta tras el sorteo. La propuesta es representarla visualmente en cada pausa del depurador — por ejemplo, como un diagrama de barras o de tarta que muestre el "peso" relativo de cada reacción candidata en ese instante, actualizándose a medida que se avanza paso a paso. Esto convierte una abstracción matemática (la distribución de propensidades) en algo que se puede ver cambiar en tiempo real, y ayuda a entender visualmente por qué una reacción "gana" a las demás en un instante dado. Coste computacional bajo: los datos ya existen, solo falta representarlos.

b) Retroceso (undo de reacciones) Permitir deshacer la última reacción disparada — o varias — y volver a un estado anterior exacto del sistema. Es la funcionalidad más compleja de las tres, por una razón de fondo: el algoritmo de Gillespie es un proceso de Markov con generación de números aleatorios; deshacer una reacción no es solo revertir el cambio de cantidades moleculares, sino también restaurar el estado del generador aleatorio y el reloj de simulación al instante exacto anterior, para que un "step" posterior sea reproducible. La forma más robusta de conseguirlo sin reescribir el motor es que el depurador mantenga snapshots periódicos del estado completo (especies + semilla del generador aleatorio + tiempo) en vez de intentar invertir matemáticamente cada reacción, y reconstruya el punto exacto navegando desde el snapshot más cercano. Es una funcionalidad de alto valor diferenciador (ningún whole-cell model existente la ofrece) pero de complejidad no trivial — se aborda después de que el modo paso a paso básico esté sólido y bien probado.

4.5 Análisis causal retrospectivo ("modo forense")

El depurador paso a paso (secciones 4.1–4.4) y el análisis causal retrospectivo responden a preguntas complementarias pero opuestas en dirección:

Depurador paso a paso:          ¿Qué va a pasar?     (hacia adelante)
Análisis causal retrospectivo:  ¿Por qué pasó?        (hacia atrás)

Dado un evento ya ocurrido en el log (por ejemplo, "la célula se quedó sin ATP en el segundo 812"), este modo permite hacer preguntas retrospectivas sobre la cadena causal que llevó a ese estado: qué reacciones consumieron más ATP en la ventana previa, qué gen dejó de expresarse justo antes, y en qué punto la trayectoria empezó a divergir del comportamiento "esperado". A diferencia del modo paso a paso (que se recorre hacia adelante y ejecuta simulación nueva), el análisis retrospectivo parte del final y reconstruye la explicación hacia atrás sobre el log ya registrado — no necesita ejecutar nada de nuevo, solo consultar y correlacionar el historial.

5. Arquitectura del sistema
Capa 1 — Motor de simulación (propio, genérico y reutilizable)
Gestor de especies: catálogo de moléculas presentes (ADN, ARNm, proteínas, metabolitos, complejos) con sus números de copias actuales.
Gestor de reacciones: cada reacción define reactivos, productos y una función de propensión (tasa dependiente del estado).
Bucle SSA propio: implementación del algoritmo de Gillespie (método directo) escrita en Python, estructurada explícitamente como generador/iterador para poder pausarse en cada paso.
Capa de depuración: envuelve el bucle SSA anterior, añadiendo modo paso a paso, evaluación de breakpoints, y registro del historial completo de eventos.
Registro de trayectorias: histórico tiempo × especie × cantidad, derivado del log de eventos, en pandas.DataFrame.
Capa 2 — Módulos biológicos

Construidos sobre la capa 1, cada uno aporta especies y reacciones al motor:

GenomeModule — Representa el ADN (secuencia + anotación de genes) usando BioPython. Genera automáticamente las reacciones de transcripción disponibles a partir de la secuencia real.
TranscriptionModule — Unión de ARN-polimerasa al promotor → elongación → terminación → liberación de ARNm. Consume NTPs.
TranslationModule — Unión de ribosoma al ARNm → elongación codón a codón (código genético real vía BioPython) → terminación → liberación de proteína. Consume aminoácidos y GTP/ATP.
DegradationModule — Degradación estocástica de ARNm y proteínas según tasas de vida media conocidas.
MetabolismModule — Red mínima de reacciones (glucólisis simplificada) que repone los precursores consumidos por transcripción/traducción. Punto de acoplamiento energético del modelo.
ReplicationModule — Dispara replicación del ADN y división celular cuando se cumplen condiciones, con segregación binomial de especies entre células hijas.
Capa 3 — Orquestador Cell + interfaz de depuración

Clase que compone los módulos activos, mantiene el estado global, y expone tanto ejecución normal como depuración:

python
cell = Cell()
cell.load_genome("mi_secuencia.fasta")

# Ejecución normal
cell.run(duration=3600)

# Ejecución depurada
debugger = cell.debug()
debugger.set_breakpoint(lambda state: state["ATP"] < 50)
debugger.step()                 # avanza una reacción
debugger.step(n=10)             # avanza 10 reacciones
debugger.run_until_breakpoint()
print(debugger.current_event)   # última reacción disparada y su contexto
print(debugger.pending_propensities())  # qué otras reacciones podían haber ocurrido
6. Genoma mínimo del MVP — y diseño para escalar
6.1 Punto de partida

Genoma de juguete de 5-8 genes esenciales, inspirado en categorías funcionales de JCVI-syn3A:

1–2 genes de proteínas ribosomales (maquinaria de traducción).
1 gen metabólico (enzima glicolítica simplificada).
1 gen regulador (represor o activador simple), para observar dinámica no trivial — y para tener un caso de uso claro del depurador (ej. "¿por qué se detuvo la transcripción de este gen?").
(Opcional) tratar la ARN-polimerasa como maquinaria preexistente, no codificada por gen, para simplificar la Fase 1.
6.2 Requisito de diseño: escalar sin reescribir

Aunque el MVP arranca con un puñado de genes, el objetivo es que el software sirva más adelante para estudios serios con genomas reales de cientos de genes (por ejemplo, aproximarse al genoma completo de JCVI-syn3A, ~473 genes), sin que eso requiera rediseñar el motor. Esto se consigue con tres decisiones concretas desde el principio:

El genoma nunca se codifica a mano en Python. Se carga siempre desde un archivo externo (FASTA para la secuencia + una anotación en formato tabular — GFF3 o un CSV simple — para posición, promotor y marco de lectura de cada gen). Añadir genes es editar el archivo de entrada, no tocar código. Esto también es lo que permite, en el futuro, importar genomas anotados reales descargados de bases de datos públicas (NCBI, EBI) sin ningún trabajo adicional.
Las reacciones se generan programáticamente a partir de la anotación, no se declaran una a una. El GenomeModule recorre la lista de genes anotados y genera automáticamente las reacciones de transcripción/traducción correspondientes a cada uno, siguiendo una plantilla común. Pasar de 8 a 400 genes no implica escribir 50 veces más código — implica que el mismo bucle de generación se ejecuta más veces.
El motor SSA propio (sección 4) no asume un número fijo de especies o reacciones. Las estructuras de datos (gestor de especies, gestor de reacciones) están indexadas dinámicamente, así que el tamaño del sistema es un parámetro de entrada, no una constante del código.
6.3 El límite real: rendimiento, no arquitectura

Como se apuntó en la sección 12, el bucle SSA propio en Python puro es más lento que una implementación optimizada en C++ como la de GillesPy2. Con 5-8 genes esto es irrelevante; con cientos de genes y miles de reacciones, el coste por paso de simulación puede volverse notable en un equipo doméstico. Este es un problema conocido y con soluciones estándar, a aplicar solo si se llega a ese punto, sin que afecten al diseño ni al código de los módulos biológicos:

Perfilar primero, optimizar después: identificar qué parte del bucle (cálculo de propensidades, muestreo aleatorio) es el cuello de botella real antes de optimizar nada.
Vectorización con numpy del cálculo de propensidades cuando el número de reacciones sea grande, en vez de iterarlas una a una en Python puro.
Tau-leaping en vez de SSA exacto para las especies muy abundantes, reduciendo el número de pasos necesarios sin perder la naturaleza estocástica del modelo (es el mismo principio que usa GillesPy2 internamente).
Extensión selectiva en Cython del núcleo del bucle SSA (no de la capa de depuración), si tras perfilar se confirma que ahí está el cuello de botella — manteniendo el resto del código en Python puro.
Reutilizar el formato SBML (ya contemplado en la sección 7) como puente: si en algún estudio serio conviene delegar la simulación en bruto a un solver externo más rápido (como el de GillesPy2) para un genoma grande, el modelo puede exportarse a SBML sin perder compatibilidad — reservando el motor propio y su depurador para el trabajo exploratorio y de inspección detallada.
7. Formatos, estándares y librerías
Necesidad	Herramienta	Propia o de terceros
Bucle de simulación SSA + depurador	Implementación propia	Propia (es el diferenciador)
Secuencias ADN/ARN, código genético, traducción	BioPython	Terceros
Redes de reacciones interoperables	SBML vía libsbml	Terceros
Cálculo numérico general	numpy	Terceros
Registro y análisis de trayectorias	pandas	Terceros
Visualización	matplotlib, plotly	Terceros
Interfaz de depuración (si se hace interactiva en terminal o notebook)	rich / ipywidgets (a evaluar)	Terceros
Metabolismo a escala genómica (futuro)	COBRApy	Terceros
7.1 Esquema de identificación de especies

Para que el modelo sea comparable e interoperable con trabajos publicados, cada especie molecular (gen, ARNm, proteína, metabolito) necesita un identificador único y consistente desde el principio — no como algo a resolver después, porque cambiarlo más tarde obligaría a tocar todos los módulos.

ID interno estable: cada especie recibe un identificador corto y único dentro del modelo (ej. ATP, mRNA_geneA, Prot_geneA), usado internamente por el motor SSA y por todos los módulos, y derivado directamente del gen de origen para que sea trazable a simple vista.
Alcance en el MVP: en esta fase el esquema cubre únicamente la identificación interna. La correspondencia con ontologías externas (BioCyc, KEGG u otras) queda fuera de alcance por ahora — se revisará si el proyecto llega a un punto en que comparar directamente con modelos publicados se vuelva necesario.

Este esquema es, además, un prerrequisito silencioso para poder exportar e importar modelos en SBML de forma correcta (sección 7), ya que SBML también exige identificadores de especie únicos.

8. Roadmap por fases
Fase	Objetivo	Criterio de éxito
0	Motor SSA propio (como generador/iterador) + clase Cell esqueleto	Simular una reacción de prueba (A + B → C) reacción a reacción, verificando que la distribución estocástica converge al comportamiento esperado
0.5	Capa de depuración sobre el motor de la Fase 0	Poder pausar, hacer step, poner un breakpoint simple, e inspeccionar el estado — sobre el sistema de prueba
1	Transcripción/traducción de 1 gen	Observar aparición y degradación de ARNm y proteína, y poder depurar por qué ocurrió cada evento
2	Varios genes + regulación simple	Dinámica de un represor afectando otro gen; usar el depurador para explicar un caso concreto (ej. "por qué se bloqueó la transcripción en el instante t")
3	Metabolismo mínimo acoplado (ATP)	La expresión génica se ralentiza si el ATP cae; breakpoint de ejemplo: pausar cuando ATP < umbral
4	Replicación de ADN + división celular	Ciclo celular completo con reparto estocástico de moléculas entre células hijas
5	Interfaz de depuración pulida + visualización básica	Depurador usable cómodamente (terminal enriquecida o notebook), gráficas de trayectorias enlazadas con el log de eventos
5.1	Visualización en tiempo real del espacio de propensidades	Ver, en cada pausa, el peso relativo de todas las reacciones candidatas, no solo la disparada
5.2	Análisis causal retrospectivo (sección 4.5)	Dado un evento del log, reconstruir hacia atrás qué reacciones y qué genes contribuyeron a ese estado
5.3	Retroceso (undo de reacciones)	Deshacer una o varias reacciones y volver a un estado anterior reproducible, vía snapshots periódicos del estado + semilla aleatoria
6	Prueba de escalabilidad	Cargar un genoma anotado real de decenas/cientos de genes (descargado de una base de datos pública) y confirmar que el GenomeModule genera las reacciones correctamente sin cambios de código; perfilar el rendimiento resultante
9. Extensiones futuras consideradas (no en el MVP)

Durante el diseño se evaluaron otros tres ángulos de innovación, aparcados deliberadamente para no dispersar el esfuerzo inicial:

Narrativa causal en lenguaje natural: traducir el log de eventos del depurador a explicaciones textuales del tipo "el represor X bloqueó la transcripción de Y, lo que redujo el ATP disponible en un Z%". Es una extensión natural de la Fase 0.5, ya que reutiliza directamente el log de eventos del depurador. Riesgo bajo si se implementa primero como reglas/plantillas, antes de considerar un LLM.
Construcción de modelos asistida por LLM: uso de un LLM para proponer automáticamente redes de reacciones y parámetros cinéticos a partir de anotaciones génicas, contrastando contra bases de datos como EcoCyc/BioCyc. Aparcado por el riesgo de que el LLM "alucine" parámetros biológicos sin trazabilidad — abordar solo con validación estricta de fuentes.
Control de versiones tipo "git" para modelos biológicos: tratar el genoma y la red de reacciones como código versionable (commits, ramas, diffs), respondiendo a una necesidad explícitamente señalada por la comunidad de whole-cell modeling. Aparcado hasta tener un modelo estable que merezca la pena versionar.
10. Estrategia de código abierto y comunidad

El proyecto se desarrollará públicamente en GitHub desde la Fase 0, no como un paso posterior de "publicación" al final. Motivos:

Trazabilidad desde el origen: cada decisión de diseño de este white paper puede enlazarse a commits e issues concretos, lo cual es coherente con la propia filosofía del depurador (todo debe ser inspeccionable y reproducible).
Validación externa: al tratarse de un modelo biofísico riguroso, exponer el código permite que otros revisen las tasas cinéticas usadas, los supuestos del modelo, y detecten errores — algo especialmente valioso dado que no formamos parte de un laboratorio con revisión por pares interna.
Responde a una necesidad ya señalada por la comunidad: como se documentó en la sección 1.2, la comunidad de whole-cell modeling ha pedido explícitamente mejores herramientas de trazabilidad y colaboración abierta para este tipo de modelos.
Reutilización futura: si más adelante se retoma la idea aparcada de "control de versiones tipo git para modelos biológicos" (sección 9), tener el proyecto ya versionado en GitHub desde el principio es una base natural sobre la que construir esa funcionalidad.

Elementos mínimos del repositorio:

Licencia abierta (a decidir — MIT o GPLv3 son las más comunes en este ecosistema; GPLv3 es la usada por E-Cell4, lo cual facilitaría eventuales colaboraciones cruzadas).
Este white paper como documento de diseño versionado (docs/ o DESIGN.md), actualizado a medida que el proyecto evolucione.
Issues abiertos para cada fase del roadmap, de forma que el progreso sea visible públicamente.
README con instrucciones de instalación y un ejemplo mínimo reproducible (el sistema de prueba A + B → C de la Fase 0), para que cualquiera pueda validar el motor SSA de forma inmediata.
Tests automatizados desde el principio (validación estadística contra soluciones analíticas, mencionada en la sección 11), configurados como integración continua (GitHub Actions) para que cada cambio se valide automáticamente.
11. Validación y criterios de éxito
Motor SSA propio: comparar contra soluciones analíticas de sistemas simples (nacimiento-muerte de una especie) para verificar que la implementación de Gillespie es correcta antes de añadir la capa de depuración.
Transcripción/traducción: usar tasas tomadas de literatura real (valores publicados para E. coli o M. genitalium), no inventadas.
Acoplamiento energético: verificar que un déficit de ATP ralentiza efectivamente la expresión génica.
Ciclo celular: verificar que el reparto de moléculas entre células hijas sigue una distribución binomial.
Depurador: verificar que el modo paso a paso produce exactamente la misma distribución estadística de resultados que el modo de ejecución normal (el depurador no debe alterar la física del sistema, solo la forma de observarlo).
12. Limitaciones conocidas
No espacial: se ignoran gradientes de concentración, difusión y geometría celular real.
Sin mutación: el ADN es estático durante toda la simulación (decisión deliberada, no limitación técnica).
Metabolismo simplificado: una vía mínima, no una red genómica completa.
Escala del genoma: 5–8 genes iniciales, lejos de los 473 de una célula mínima real; escalar es principalmente un problema de rendimiento del bucle SSA propio, a resolver con optimización (ej. Cython) si se llega a ese punto.
Motor propio vs. librerías optimizadas: al escribir el bucle SSA nosotros mismos en Python puro, será más lento que soluciones en C++ como GillesPy2 para modelos grandes. Aceptable para el MVP (pocas reacciones); revisar si se escala mucho el genoma.
13. Referencias clave
Karr, J.R. et al. (2012). A Whole-Cell Computational Model Predicts Phenotype from Genotype. Cell, 150(2), 389–401.
Thornburg, Z.R. et al. (2022). Fundamental behaviors emerge from simulations of a living minimal cell. Cell, 185(2), 345–360.
Gillespie, D.T. (1977). Exact stochastic simulation of coupled chemical reactions. The Journal of Physical Chemistry, 81(25), 2340–2361.
Hucka, M. et al. (2003). The systems biology markup language (SBML). Bioinformatics, 19(4), 524–531.
Stumpf, M.P.H. (2021). Statistical and computational challenges for whole cell modelling. Current Opinion in Systems Biology, 26, 58–63.
Karr, J.R. et al. (2017). A blueprint for human whole-cell modeling. (Propuesta de estándares y control de versiones para la comunidad).
Abel, J.H. et al. / GillesPy2 documentation (StochSS project) — referencia de por qué las librerías SSA existentes no exponen ejecución paso a paso de forma nativa.
