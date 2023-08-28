FROM python:3

RUN groupadd -r snake && useradd -r -g snake snake

WORKDIR /home/snake

COPY Pipfile* ./
RUN pip install --disable-pip-version-check --no-cache-dir pipenv && pipenv install --system --deploy

USER snake

CMD [ "python", "./gen_img_from_lay.py" ]
