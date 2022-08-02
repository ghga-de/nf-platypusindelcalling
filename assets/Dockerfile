FROM nfcore/base:1.14
LABEL authors="Kuebra narci" \
      description="Docker image containing all software requirements for the IndelCalling pipeline RunTinda"

# Install the conda environment
COPY tindaperlenvironment.yml /
RUN conda env create --quiet -f /tindaperlenvironment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/INDELCALL/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name INDELCALL > indelcall_runtinda_v1.yml
