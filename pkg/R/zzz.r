.onAttach <- function(libname, pkgname){
  pkg.version <- packageDescription("ggene", fields = "Version")
  startup.txt <- paste("\nggene version ", pkg.version, " is loaded\n", sep="")
  packageStartupMessage(startup.txt)
}
