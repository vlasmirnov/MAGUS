FROM python:3.10.2-slim
WORKDIR /usr/src/app
COPY . .
RUN apt-get update
RUN apt-get install libgomp1
RUN pip install dendropy
RUN find tools -maxdepth 5 -type f ! -name "*.*" ! -name "README" | xargs -I "{}" chmod +x {}
ENTRYPOINT [ "python3", "magus.py" ]