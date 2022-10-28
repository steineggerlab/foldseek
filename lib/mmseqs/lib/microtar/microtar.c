/*
 * Copyright (c) 2017 rxi
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>

#include "microtar.h"

typedef struct {
  char name[100];
  char mode[8];
  char owner[8];
  char group[8];
  char size[12];
  char mtime[12];
  char checksum[8];
  char type;
  char linkname[100];
  char _padding[255];
} mtar_raw_header_t;


static unsigned round_up(unsigned n, unsigned incr) {
  return n + (incr - n % incr) % incr;
}

static unsigned checksum(const mtar_raw_header_t* rh) {
  unsigned i;
  unsigned char *p = (unsigned char*) rh;
  unsigned res = 256;
  for (i = 0; i < offsetof(mtar_raw_header_t, checksum); i++) {
    res += p[i];
  }
  for (i = offsetof(mtar_raw_header_t, type); i < sizeof(*rh); i++) {
    res += p[i];
  }
  return res;
}

const char* mtar_strerror(int err) {
  switch (err) {
    case MTAR_ESUCCESS     : return "success";
    case MTAR_EFAILURE     : return "failure";
    case MTAR_EOPENFAIL    : return "could not open";
    case MTAR_EREADFAIL    : return "could not read";
    case MTAR_ESEEKFAIL    : return "could not seek";
    case MTAR_EBADCHKSUM   : return "bad checksum";
    case MTAR_ENULLRECORD  : return "null record";
    case MTAR_ENOTFOUND    : return "file not found";
  }
  return "unknown error";
}

static int file_read(mtar_t *tar, void *data, size_t size) {
  size_t res = fread(data, 1, size, (FILE*)tar->stream);
  return (res == size) ? MTAR_ESUCCESS : MTAR_EREADFAIL;
}

static int file_seek(mtar_t *tar, long offset, int whence) {
  int res = fseek((FILE*)tar->stream, offset, whence);
  return (res == 0) ? MTAR_ESUCCESS : MTAR_ESEEKFAIL;
}

static int file_close(mtar_t *tar) {
  fclose((FILE*)tar->stream);
  return MTAR_ESUCCESS;
}

int mtar_open(mtar_t *tar, const char *filename) {
  /* Init tar struct and functions */
  memset(tar, 0, sizeof(*tar));
  tar->read = file_read;
  tar->seek = file_seek;
  tar->close = file_close;

  /* Open file */
  tar->stream = fopen(filename, "rb");
  if (!tar->stream) {
    return MTAR_EOPENFAIL;
  }

  /* Return ok */
  return MTAR_ESUCCESS;
}

int mtar_close(mtar_t *tar) {
  return tar->close(tar);
}

int mtar_read_header(mtar_t *tar, mtar_header_t *h) {
  mtar_raw_header_t rh;
  /* Read raw header */
  int err = tar->read(tar, &rh, sizeof(rh));
  if (err) {
    return err;
  }

  /* Load raw header into header struct and return */
  unsigned chksum1, chksum2;

  /* If the checksum starts with a null byte we assume the record is NULL */
  if (*(rh.checksum) == '\0') {
      return MTAR_ENULLRECORD;
  }

  /* Build and compare checksum */
  chksum1 = checksum(&rh);
  chksum2 = strtoul(rh.checksum, NULL, 8);
  if (chksum1 != chksum2) {
      return MTAR_EBADCHKSUM;
  }

  /* Load raw header into header */
  h->mode =  strtoul(rh.mode,  NULL, 8);
  h->owner = strtoul(rh.owner, NULL, 8);
  h->size =  strtoul(rh.size,  NULL, 8);
  h->mtime = strtoul(rh.mtime, NULL, 8);
  h->type = rh.type;
  if (h->type != MTAR_TGNU_LONGNAME && h->type != MTAR_TGNU_LONGLINK) {
    memcpy(h->name, rh.name, sizeof(h->name));
    h->name[sizeof(h->name)-1] = '\0';
    memcpy(h->linkname, rh.linkname, sizeof(h->linkname));
    h->linkname[sizeof(h->linkname)-1] = '\0';
  } else {
    h->name[0] = '\0';
    h->linkname[0] = '\0';
  }
  tar->curr_size = h->size;

  return MTAR_ESUCCESS;
}

int mtar_skip_data(mtar_t *tar) {
  int n = round_up(tar->curr_size, 512);
  int err = tar->seek(tar, n, SEEK_CUR);
  if (err) {
      return err;
  }
  return MTAR_ESUCCESS;
}

int mtar_read_data(mtar_t *tar, void *ptr, size_t size) {
  /* Read data */
  int err = tar->read(tar, ptr, size);
  if (err) {
    return err;
  }
  int n = round_up(tar->curr_size, 512) - tar->curr_size;
  err = tar->seek(tar, n, SEEK_CUR);
  if (err) {
      return err;
  }
  return MTAR_ESUCCESS;
}
