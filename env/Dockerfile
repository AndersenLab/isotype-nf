FROM python:3-slim
# RUN apk --no-cache --update-cache add gcc gfortran python python-dev py-pip build-base wget freetype-dev libpng-dev openblas-dev
# RUN ln -s /usr/include/locale.h /usr/include/xlocale.h
RUN pip install --no-cache-dir numpy scipy matplotlib cairocffi
RUN apt-get --allow-releaseinfo-change update && \
   apt-get install -y procps libgtk-3-dev wget && \
	rm -rf /var/lib/apt/lists/*