.onLoad <- function(libname, pkgname) {
    utils::globalVariables("org.Sco.egGO")
    utils::globalVariables("ppiPreEnv")
    utils::globalVariables("org.Sco.egGO2EG")

   	.initial()
}
