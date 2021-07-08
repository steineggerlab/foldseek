FROM debian:stable-slim as foldseek-builder

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
    build-essential cmake xxd git zlib1g-dev libbz2-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/foldseek
ADD . .

RUN mkdir -p build_sse/bin && mkdir -p build_avx/bin

WORKDIR /opt/foldseek/build_sse
RUN cmake -DHAVE_SSE4_1=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    make -j $(nproc --all) && make install;

WORKDIR /opt/foldseek/build_avx
RUN cmake -DHAVE_AVX2=1 -DHAVE_MPI=0 -DHAVE_TESTS=0 -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=. ..; \
    make -j $(nproc --all) && make install;

FROM debian:stable-slim
MAINTAINER Martin Steinegger <martin.steinegger@snu.ac.kr>

RUN apt-get update && apt-get upgrade -y && apt-get install -y \
     gawk bash grep wget libstdc++6 libgomp1 zlib1g libbz2-1.0 \
    && rm -rf /var/lib/apt/lists/*

COPY --from=foldseek-builder /opt/foldseek/build_sse/bin/foldseek /usr/local/bin/foldseek_sse42
COPY --from=foldseek-builder /opt/foldseek/build_avx/bin/foldseek /usr/local/bin/foldseek_avx2
RUN echo '#!/bin/bash\n\
if $(grep -q -E "^flags.+avx2" /proc/cpuinfo); then\n\
    exec /usr/local/bin/foldseek_avx2 "$@"\n\
else\n\
    exec /usr/local/bin/foldseek_sse42 "$@"\n\
fi' > /usr/local/bin/foldseek
RUN chmod +x /usr/local/bin/foldseek

ENTRYPOINT ["/usr/local/bin/foldseek"]
