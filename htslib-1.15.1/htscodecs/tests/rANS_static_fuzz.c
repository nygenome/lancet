/* Fuzz testing target. */
/*
 * Copyright (c) 2019,2020 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*
Local instructions: compile, from a build subdir, with
/software/badger/opt/llvm/7.0.0/bin/clang -O3 -g ../../tests/rANS_static_fuzz.c -I../.. ../../htscodecs/rANS_static.c  -pthread -fsanitize=fuzzer,address /software/badger/opt/gcc/8.1.0/lib64/libstdc++.a -DFUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION

(This bizarrity is because our local clang install wasn't built with
C++ support.)

Run with:
    export ASAN_OPTIONS=allow_addr2line=true
    ./a.out -rss_limit_mb=8000 corpus
or 
    ./a.out -rss_limit_mb=8000 -detect_leaks=0 corpus

I generated corpus as a whole bunch of precompressed tiny inputs from
tests/dat/q4 for different compression modes.

For debugging purposes, we can compile a non-fuzzer non-ASAN build using
-DNOFUZZ which creates a binary we can debug on any libfuzzer generated
output using valgrind.  (The rans4x8 command line test won't quite work as
it's a slightly different input format with explicit sizes in the binary
stream.)
*/

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <assert.h>
#include <string.h>
#include <sys/time.h>

#include "htscodecs/rANS_static.h"
#include "htscodecs/rANS_static.c"

int LLVMFuzzerTestOneInput(uint8_t *in, size_t in_size) {
    unsigned int uncomp_size;
    unsigned char *uncomp = rans_uncompress(in, in_size, &uncomp_size);
    if (uncomp)
	free(uncomp);
    
    return 0;
}

#ifdef NOFUZZ
#include <stdio.h>
#include <fcntl.h>
#include <errno.h>

#define BS 1024*1024
static unsigned char *load(char *fn, uint64_t *lenp) {
    unsigned char *data = NULL;
    uint64_t dsize = 0;
    uint64_t dcurr = 0;
    signed int len;
    int fd = open(fn, O_RDONLY);

    do {
	if (dsize - dcurr < BS) {
	    dsize = dsize ? dsize * 2 : BS;
	    data = realloc(data, dsize);
	}

	len = read(fd, data + dcurr, BS);
	if (len > 0)
	    dcurr += len;
    } while (len > 0);

    if (len == -1) {
	perror("read");
    }

    close(fd);
    *lenp = dcurr;
    return data;
}

int main(int argc, char **argv) {
    uint64_t in_size;
    unsigned char *in = load(argv[1], &in_size);
    unsigned int uncomp_size;
    unsigned char *uncomp = rans_uncompress(in, in_size, &uncomp_size);
    if (uncomp)
	free(uncomp);

    free(in);
    
    return 0;
}
#endif
