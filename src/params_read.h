// JSON parameter reader for BioReactor.c.
// Wraps jsmn (MIT, Serge Zaitsev) to parse params.json into a flat C struct.
// All physics parameters live here; BioReactor.c reads only from BioreactorParams.
//
// Canonical params.json schema (see design spec):
//   fidelity, omega_b, n_harmonics,
//   theta_max[3], phi_angular[3],
//   omega_h, amplitude_h[3], phi_horizontal[3],
//   geometry.{a, b, n}, fill_level
//
// phi_angular[0] is forced to 0.0 at read time (time-origin reference, not a free parameter).

#ifndef PARAMS_READ_H
#define PARAMS_READ_H

#define JSMN_STATIC
#include "jsmn.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N_MAX 3  // maximum number of harmonics; vectors always padded to this length

typedef struct {
  int    fidelity;
  double omega_b;
  int    n_harmonics;
  double theta_max[N_MAX];
  double phi_angular[N_MAX];    // [0] always 0
  double omega_h;
  double amplitude_h[N_MAX];
  double phi_horizontal[N_MAX];
  double geometry_a;
  double geometry_b;
  double geometry_n;
  double fill_level;
  double t_end;                 // simulation end time [non-dim]; default 250.0
  int    n_mix_cycles;          // rocking cycles before oxygen/tracer start; default 80
} BioreactorParams;

// ── helpers ──────────────────────────────────────────────────────────────────

static int jsoneq(const char *json, jsmntok_t *tok, const char *s) {
  return tok->type == JSMN_STRING
      && (int)strlen(s) == tok->end - tok->start
      && strncmp(json + tok->start, s, tok->end - tok->start) == 0;
}

static double tok_double(const char *json, jsmntok_t *t) {
  char buf[64];
  int len = t->end - t->start;
  if (len >= (int)sizeof(buf)) len = (int)sizeof(buf) - 1;
  memcpy(buf, json + t->start, len);
  buf[len] = '\0';
  return atof(buf);
}

static int tok_int(const char *json, jsmntok_t *t) {
  return (int)tok_double(json, t);
}

// Read a JSON array of up to N_MAX doubles; remaining entries stay at 0.
static void tok_array(const char *json, jsmntok_t *tokens, int arr_idx,
                      double *out, int n_max) {
  jsmntok_t *arr = &tokens[arr_idx];
  int count = arr->size < n_max ? arr->size : n_max;
  for (int i = 0; i < count; i++)
    out[i] = tok_double(json, &tokens[arr_idx + 1 + i]);
}

// ── public API ───────────────────────────────────────────────────────────────

static BioreactorParams params_read(const char *path) {
  BioreactorParams p = {0};  // zero-initialise; pads harmonic vectors to 0
  p.n_harmonics   = 1;         // default: single sinusoid (pure rocking)
  p.t_end         = 250.0;    // default if not present in params.json
  p.n_mix_cycles  = 80;       // default: 80 rocking cycles (upstream hardcoded value)

  FILE *fp = fopen(path, "r");
  if (!fp) {
    fprintf(stderr, "params_read: cannot open '%s'\n", path);
    exit(1);
  }
  fseek(fp, 0, SEEK_END);
  long sz = ftell(fp);
  rewind(fp);
  char *json = (char *)malloc(sz + 1);
  fread(json, 1, sz, fp);
  fclose(fp);
  json[sz] = '\0';

  jsmn_parser parser;
  jsmn_init(&parser);
  // count tokens first
  int ntok = jsmn_parse(&parser, json, sz, NULL, 0);
  if (ntok < 0) {
    fprintf(stderr, "params_read: jsmn_parse failed (%d) on '%s'\n", ntok, path);
    exit(1);
  }
  jsmntok_t *tokens = (jsmntok_t *)malloc(ntok * sizeof(jsmntok_t));
  jsmn_init(&parser);
  jsmn_parse(&parser, json, sz, tokens, ntok);

  // top-level object
  for (int i = 1; i < ntok; i++) {
    if (jsoneq(json, &tokens[i], "fidelity"))
      p.fidelity = tok_int(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "omega_b"))
      p.omega_b = tok_double(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "n_harmonics"))
      p.n_harmonics = tok_int(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "omega_h"))
      p.omega_h = tok_double(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "fill_level"))
      p.fill_level = tok_double(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "t_end"))
      p.t_end = tok_double(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "n_mix_cycles"))
      p.n_mix_cycles = tok_int(json, &tokens[++i]);
    else if (jsoneq(json, &tokens[i], "theta_max")) {
      tok_array(json, tokens, ++i, p.theta_max, N_MAX);
      i += tokens[i].size;
    }
    else if (jsoneq(json, &tokens[i], "phi_angular")) {
      tok_array(json, tokens, ++i, p.phi_angular, N_MAX);
      i += tokens[i].size;
    }
    else if (jsoneq(json, &tokens[i], "amplitude_h")) {
      tok_array(json, tokens, ++i, p.amplitude_h, N_MAX);
      i += tokens[i].size;
    }
    else if (jsoneq(json, &tokens[i], "phi_horizontal")) {
      tok_array(json, tokens, ++i, p.phi_horizontal, N_MAX);
      i += tokens[i].size;
    }
    else if (jsoneq(json, &tokens[i], "geometry")) {
      // geometry is a nested object; walk its key-value pairs
      i++;  // move to the object token
      int geo_size = tokens[i].size;
      for (int g = 0; g < geo_size; g++) {
        i++;
        if (jsoneq(json, &tokens[i], "a"))
          p.geometry_a = tok_double(json, &tokens[++i]);
        else if (jsoneq(json, &tokens[i], "b"))
          p.geometry_b = tok_double(json, &tokens[++i]);
        else if (jsoneq(json, &tokens[i], "n"))
          p.geometry_n = tok_double(json, &tokens[++i]);
        else
          i++;  // skip unknown geometry keys
      }
    }
  }

  // phi_angular[0] is the time-origin reference — always 0, not a free parameter
  p.phi_angular[0] = 0.0;

  free(tokens);
  free(json);
  return p;
}

#endif // PARAMS_READ_H
