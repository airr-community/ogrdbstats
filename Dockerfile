FROM immcantation/suite:4.4.0
COPY --chmod=777 run_ogrdbstats /usr/local/bin
COPY --chmod=777 run_ogrdbstats_tests /usr/local/bin
COPY --chmod=777 make_sample_data /usr/local/bin
COPY inst/extdata /data/ogrdb_extdata
COPY *.R /usr/local/bin
WORKDIR /usr/local/bin
RUN Rscript install_ogrdbstats.R


