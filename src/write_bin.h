/*!
 *  \file write_bin.h
 *  \brief      Utility functions for writing data to files
 *  \author    Damian Miralles
 *  \version   4.1a
 *  \date      Jan 23, 2018
 */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

int write_file_fl64(const char *file_path, double *data) {

  FILE *fp = fopen(file_path, "wb");
  fwrite(data, sizeof *data, 37000, fp);
  fclose(fp);

  return EXIT_SUCCESS;
}
