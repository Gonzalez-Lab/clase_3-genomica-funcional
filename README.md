Clase 3 — Genómica funcional: análisis de ATAC-seq y ChIP-seq

Curso: Descifrando el código epigenético: de la genómica funcional al análisis bioinformático de datos ómicos

Institución: ITBA · 2026

Material práctico: Gonzalez Lab

Objetivo
Esta clase introduce el análisis bioinformático de regiones genómicas enriquecidas (peaks) obtenidas mediante ATAC-seq y ChIP-seq. A partir de archivos de picos, alineamientos y resultados de análisis diferencial, se trabaja sobre el control de calidad, la anotación genómica y la interpretación de cambios entre condiciones experimentales.

En el caso de ATAC-seq, los picos representan regiones de cromatina accesible. En el caso de ChIP-seq, representan regiones enriquecidas para la proteína o modificación de histonas investigada; en este práctico se trabaja con H3K27ac.
Contenido del repositorio
Archivo	Función
[`Peaks files processing (ATAC or ChIP).R`](Peaks files processing (ATAC or ChIP\).R)	Importación de archivos narrowPeak, construcción de regiones consenso entre réplicas, anotación con ChIPseeker, distribución genómica, distancia al TSS y asociaciones pico–gen.
[`FRiP_TSS_ATAC.R`](FRiP_TSS_ATAC.R)	Evaluación de calidad de ATAC-seq: fracción de fragmentos en picos (FRiP) y enriquecimiento de señal alrededor de los sitios de inicio de transcripción (TSS).
[`DiffBind_analysis.R`](DiffBind_analysis.R)	Análisis diferencial de enriquecimiento entre las condiciones SAHA y Control, visualizaciones, anotación de regiones diferenciales y selección de picos asociados a promotores.
data.rar	Archivo comprimido con datos de trabajo. Descomprimir antes de ejecutar los scripts y comprobar las rutas de los archivos.
db_ChIP.rds	Objeto R serializado para retomar el análisis de ChIP-seq sin repetir las etapas computacionalmente más costosas, cuando sea compatible con el script y la versión instalada de los paquetes.


Flujo de trabajo
1. Procesamiento de picos: importar los narrowPeak de las réplicas, identificar regiones con soporte en al menos dos réplicas y anotar su ubicación respecto de genes y TSS.
2. Control de calidad de ATAC-seq: estimar FRiP y examinar el enriquecimiento alrededor de los TSS para cada muestra.
3. Análisis diferencial: utilizar DiffBind para cuantificar la señal sobre un conjunto común de regiones, comparar SAHA frente a Control y recuperar los picos con diferencias estadísticamente significativas.
4. Interpretación: distinguir regiones con mayor o menor señal, explorar su distribución genómica y obtener asociaciones entre picos de promotores y genes cercanos.
Importante: el criterio de soporte en 2 de 3 réplicas permite construir un conjunto de picos consenso, pero no reemplaza un análisis formal de reproducibilidad como IDR. Asimismo, asignar un pico a un gen no demuestra un cambio en su expresión.

Gonzalez Lab · ITBA · 2026
