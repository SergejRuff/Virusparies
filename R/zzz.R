.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "*************************************************************\n",
    "*   __     ___                                 _            *\n",
    "*   \\ \\   / (_)_ __ _   _ ___ _ __   __ _ _ __(_) ___  ___  *\n",
    "*    \\ \\ / /| | '__| | | / __| '_ \\ / _` | '__| |/ _ \\/ __| *\n",
    "*     \\ V / | | |  | |_| \\__ \\ |_) | (_| | |  | |  __/\\__ \\ *\n",
    "*      \\_/  |_|_|   \\__,_|___/ .__/ \\__,_|_|  |_|\\___||___/ *\n",
    "*                            |_|                            *\n",
    "*                                                           *\n",
    "*                                                           *\n",
    "*************************************************************\n",
    "*                                                           *\n",
    "*               Welcome to Virusparies!                     *\n",
    "*                                                           *\n",
    "*************************************************************\n",
    "\n",
    "Create stunning visuals for VirusHunterGatherer hittables.\n",
    "Enjoy exploring your data!\n",
    "\n"
  )

  packageStartupMessage("This is version ", utils::packageVersion(pkgname),
                        " of ", pkgname,".")
}

# Declare global variables
utils::globalVariables(c("ICTV_env"))
