FROM debian:10-slim

RUN apt update && apt install -y build-essential gcc g++ gfortran git sudo tar curl zip unzip \
        ninja-build pkg-config gdb autoconf libtool python3 \
        ca-certificates libasound2 libatk1.0-0  \
        libgstreamer-plugins-base1.0-0 libgstreamer1.0-0 libgtk-3-0 libidl-2-0 \ 
        libxft2 libxi6 libxinerama1 libxrandr2 libxrender1 libxss1 libxt6 \
        libnss3 libpam0g libpango-1.0-0 libpangocairo-1.0-0 libpangoft2-1.0-0 \
        libpython2.7 libpython3.6 libpython3.7 libselinux1 libsm6 libsndfile1 libx11-6 \
        libx11-xcb1 libxcb1 libxcomposite1 libxcursor1 libxdamage1 libxext6 libxfixes3 \
        libfontconfig1 libgconf-2-4 libgdk-pixbuf2.0-0 \
        libc6 libcairo-gobject2 libcairo2 libcups2 libdbus-1-3 libdbus-glib-1-2 \
        libxtst6 libxxf86vm1 procps sudo unzip wget xkb-data zlib1g libavcodec-extra58 libavformat58

WORKDIR /tmp

RUN wget https://www.mathworks.com/mpm/glnxa64/mpm \
        && chmod +x mpm \
        && ./mpm install  --release R2019b --products MATLAB \
        && update-alternatives --install /usr/bin/matlab matlab /usr/share/matlab/bin/matlab 10

ADD license.lic.gpg /usr/share/matlab/licenses/license.lic.gpg

ADD reinstall-cmake.sh /tmp/reinstall-cmake.sh

RUN chmod +x /tmp/reinstall-cmake.sh \
    && /tmp/reinstall-cmake.sh 3.29.6

ARG USERNAME=vscode

RUN adduser --shell /bin/bash --disabled-password --gecos "" ${USERNAME} \
    && echo "$USERNAME ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/${USERNAME} \
    && chmod 0440 /etc/sudoers.d/${USERNAME} \
    && chown -R ${USERNAME} /usr/share/matlab/licenses


ENV LD_LIBRARY_PATH=.:/usr/share/matlab/runtime/glnxa64:/usr/share/matlab/bin/glnxa64:/usr/share/matlab/sys/os/glnxa64:/usr/share/matlab/sys/opengl/lib/glnxa64
ENV LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6:/lib/x86_64-linux-gnu/libc.so.6