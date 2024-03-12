FROM python:3.11.5-bullseye

ENV PYTHONUNBUFFERED=1
ENV HOME=/code

COPY requirements.txt /tmp/
RUN pip install -r /tmp/requirements.txt

WORKDIR ${HOME}
COPY src/ ./
