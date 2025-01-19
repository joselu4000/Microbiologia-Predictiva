Autor: José Luis López Carmona
Contacto: joselu4000@gmail.com

Este repositorio responde a los códigos usados para el Trabajo Final del Máster Universitario de Matemáticas por la Universidad de Sevilla.

A continuación, procedo a describir la notación de los métodos usados para su correcta identificación (según el trabajo escrito). Refiriendo a las distintas carpetas:

**datos**:
- datos.xlsx, son los datos 1 del trabajo.
- Datos2.xlsx, son los datos 2. Atención porque contiene más de los usados.
- ComBaseExport.xlsx, los datos usados para los ejemplos de latencia.
- Los usados para el Modelo Secundario no están incluidos puesto que la autoría y propiedad son de un artículo publicado (Heredia et al., 2024). Si se quieren, para un uso científico y/o académico, pueden ponerse en contacto conmigo.

**Trifase Primario:**
- Gauss-Newton: TrifaseGNRF1True.R y TrifaseGNRFTrue.R, para datos 1 y datos 2.
- Adam: TrifaseAdamD1.R y TrifaseAdamD2.R
- ode-Adam: TrifaseANAdamD1.R y TrifaseANAdamD2.R
- ode-L-BFGS-B: TrifaseAND1.R y TrifaseANR.R
Incluido un shyni para manipulación de parámetros directa.

**Logistico Primario**:
- Gauss-Newton: GompertzGNRD1.R y GompertzGNR.R
- Adam: GompertzAdamD1.R y GompertzAdamD2.R
- ode-Adam: GompertzANAdamD1.R y GompertzANAdamD2.R
- ode-L-BFGS-B: GompertzANRD1.R y GompertzANR.R
Incluido un shyni para manipulación de parámetros directa.

**Baranyi Primario**:
- Gauss-Newton: BaranyiGND1.R y BaranyiGNR.R
- Adam: BaranyiAdamD1.R y BaranyiAdamD2.R
- ode-Adam: BaranyiANAdamD1.R y BaranyiANAdamD2.R
- ode-L-BFGS-B: BaranyiANRD1.R y BaranyiANR.R
Incluido un shyni para manipulación de parámetros directa.

**LagExpo Primario**:
- Gauss-Newton: LagExpoGND1.R y LagExpoGNR.R
- Adam: LagExpoAdamD1.R y LagExpoAdamD2.R
Incluido un shyni para manipulación de parámetros directa.

**Sensibilidad**:
- SensibilidadTRifase.R
- SensibilidadLogistico.R
- SensibilidadBaranyi.R

**Figures**, donde quedarán contenidas la mayor parte de figuras tras la carga de los códigos.

**Archivos sin carpeta**:
- figures.R, un generador de figuras.
- Secundario.R, el código usado para el Modelo Secundario

Todos los códigos están mínimamente comentados. En caso de dudas pueden ponerse en contacto conmigo.






