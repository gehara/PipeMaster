#ifndef MSABC_CAPTURE_H
#define MSABC_CAPTURE_H

#include <stdio.h>

/* Output capture for R package mode */
FILE *msABC_get_output_stream(void);
char *msABC_get_output_buffer(size_t *len);
void msABC_close_output_stream(void);
void msABC_init_output_stream(void);

/* Global state management */
extern void msABC_reset_streec_statics(void);
extern void msABC_set_jmpbuf_active(int active);
extern int msABC_main(int argc, char **argv);

/* Seed management */
extern void seedit_r(unsigned short *seedv);
extern void get_seed_r(unsigned short *seedv);

#endif /* MSABC_CAPTURE_H */
