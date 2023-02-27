FROM python:3.8

ADD main.py .

RUN pip install numpy
RUN pip install matplotlib

CMD [ "python", "./main.py" ]
