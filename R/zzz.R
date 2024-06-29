.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "*********************************************************\n",
    "*           __           _     ______                   *\n",
    "*          \\   \\       /  /   |  __  \\                  *\n",
    "*           \\   \\     /  /    | |  | |                  *\n",
    "*            \\   \\   /  /     | |__| |                  *\n",
    "*             \\   \\_/  /      |   _ /                   *\n",
    "*              \\      /       |  |                      *\n",
    "*               \\    /        |  |                      *\n",
    "*                \\__/         |__|                      *\n",
    "*                                                       * \n",
    "*********************************************************\n",
    "*                                                       *\n",
    "*               Welcome to Virusparies!                 *\n",
    "*                                                       *\n",
    "*********************************************************\n",
    "\n",
    "Create stunning visuals for VirusHunterGatherer hittables.\n",
    "Enjoy exploring your data!\n",
    "\n"
  )

  packageStartupMessage("This is version ", utils::packageVersion(pkgname),
                        " of ", pkgname,".")
}

