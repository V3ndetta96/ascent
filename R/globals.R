# ==============================================================================
# DECLARACIÓN DE VARIABLES GLOBALES
# ==============================================================================
# Esto evita las advertencias de 'R CMD check' por el uso de evaluación
# no estándar en ggplot2 (funciones aes()).

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    # Variables usadas en plot.ascfcd y plot.ascfcd_pw
    "PC1",
    "PC2",
    "State",
    "Community",
    "Leverage",
    "Species",

    # Variables usadas en los gráficos de transiciones / Delta CWM
    "Trait",
    "Delta_CWM_Std",
    "Abs_Shift"
  ))
}
