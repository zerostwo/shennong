.onAttach <- function(libname, pkgname) {
  if (!isTRUE(getOption("shennong.startup_message", TRUE))) {
    return(invisible())
  }
  if (!interactive() && !isTRUE(getOption("shennong.startup_message_always", FALSE))) {
    return(invisible())
  }

  packageStartupMessage(
    paste(
      "  ____  _                              ",
      " / ___|| |__   ___ _ __  _ __   ___  _ __   __ _",
      " \\___ \\| '_ \\ / _ \\ '_ \\| '_ \\ / _ \\| '_ \\ / _` |",
      "  ___) | | | |  __/ | | | | | | (_) | | | | (_| |",
      " |____/|_| |_|\\___|_| |_|_| |_|\\___/|_| |_|\\__, |",
      "                                           |___/ ",
      "",
      "Shennong is ready.",
      "Try `sn_list_results()`, `sn_list_palettes()`, or `sn_initialize_codex_project()`.",
      sep = "\n"
    )
  )

  invisible()
}
