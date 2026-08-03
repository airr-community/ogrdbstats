FROM immcantation/suite:4.7.0
COPY --chmod=777 run_ogrdbstats /usr/local/bin
COPY --chmod=777 run_ogrdbstats_tests /usr/local/bin
COPY --chmod=777 make_sample_data /usr/local/bin
COPY inst/extdata /data/ogrdb_extdata
COPY *.R /usr/local/bin
WORKDIR /usr/local/bin
RUN Rscript install_ogrdbstats.R
RUN curl -sSL "https://github.com/peak/s5cmd/releases/download/v2.3.0/s5cmd_2.3.0_Linux-64bit.tar.gz" | tar -xz -C /usr/local/bin s5cmd
RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip" \
    && unzip awscliv2.zip \
    && ./aws/install \
    && rm -rf awscliv2.zip aws

