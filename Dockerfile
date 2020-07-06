# Build modwt
FROM alpine:3.7 as modwt-build
RUN apk add --no-cache \
  g++ \
  git \
  make
RUN git clone https://github.com/StamLab/modwt.git \
      && cd modwt \
      && git checkout 28e9f479c737836ffc870199f2468e30659ab38d \
      && make

# Build samtools
FROM alpine:3.7 as samtools-build
RUN apk add --no-cache \
    autoconf \
    g++ \
    git \
    bzip2-dev \
    make \
    ncurses-dev \
    xz-dev \
    zlib-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.7/samtools-1.7.tar.bz2 \
      && tar xf samtools-1.7.tar.bz2 \
      && cd samtools-1.7 \
      && make install

# Build hotspot2
FROM alpine:3.7 as hotspot2-build
RUN apk add --no-cache \
      bash \
      g++ \
      make
WORKDIR /hotspot2
ADD . .
RUN make

# Build bedGraphToBigWig
FROM alpine:3.7 as kentutils-build
RUN apk add --no-cache \
      g++ \
      gcc \
      git \
      libpng-dev \
      make \
      mysql-dev \
      zlib-dev
# Get an archive of kentUtils, remove a file that doesn't build, and compile
RUN wget https://github.com/ENCODE-DCC/kentUtils/archive/v302.0.0.tar.gz \
      && tar xf v302.0.0.tar.gz \
      && cd kentUtils-302.0.0 \
      && sed -i 's/fof.o //' src/lib/makefile \
      && make

# Build the final container
FROM alpine:3.7 as hotspot2
# Install dynamic libraries
RUN apk add --no-cache \
      bash \
      bc \
      bzip2-dev \
      coreutils \
      ncurses \
      xz-dev \
      zlib-dev
# Get bedops
RUN wget -O - https://github.com/bedops/bedops/releases/download/v2.4.31/bedops_linux_x86_64-v2.4.31.tar.bz2 \
  | tar -C /usr/local -xjf -
# Get built files
COPY --from=kentutils-build /kentUtils-302.0.0/bin/bedGraphToBigWig /usr/local/bin/
COPY --from=kentutils-build /kentUtils-302.0.0/bin/bedToBigBed /usr/local/bin/

COPY --from=modwt-build /modwt/bin/ /usr/local/bin/

COPY --from=samtools-build /usr/local/bin /usr/local/bin/

COPY --from=hotspot2-build /hotspot2/bin/ /usr/local/bin/
COPY --from=hotspot2-build /hotspot2/scripts/ /usr/local/bin/
