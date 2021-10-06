void R_init_cucumber_librust(void *dll);
void R_init_cucumber(void *dll) { R_init_cucumber_librust(dll); }
