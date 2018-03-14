# Build modwt
FROM alpine as modwt-build
RUN apk update
RUN apk add \
  g++ \
  git \
  make
RUN git clone https://github.com/StamLab/modwt.git \
      && cd modwt \
      && git checkout 28e9f479c737836ffc870199f2468e30659ab38d \
      && make

# Build samtools
FROM alpine as samtools-build
RUN apk update \
    && apk add \
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
FROM alpine as hotspot2-build
RUN apk update
RUN apk add \
      bash \
      g++ \
      make
WORKDIR /hotspot2
ADD . .
RUN make

# Build the final container
FROM alpine as hotspot2
# Install dynamic libraries
RUN apk update \
      && apk add \
      bash \
      bzip2-dev \
      ncurses \
      xz-dev \
      zlib-dev
# Get bedops
RUN wget -O - https://github.com/bedops/bedops/releases/download/v2.4.31/bedops_linux_x86_64-v2.4.31.tar.bz2 \
  | tar -C /usr/local -xjf -
# Get bedgraph2BigWig
RUN cd /usr/local/bin && wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
# Get built files
COPY --from=modwt-build /modwt/bin/ /usr/local/bin
COPY --from=samtools-build /usr/local/bin /usr/local/bin
COPY --from=hotspot2-build /hotspot2/bin/ /usr/local/bin/
COPY --from=hotspot2-build /hotspot2/scripts/ /usr/local/bin/
