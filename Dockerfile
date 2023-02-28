FROM python:3.8

ADD main.py .
ADD argons.xyz .

RUN pip install numpy
RUN pip install matplotlib
RUN pip install xyz-py
RUN pip install scipy

CMD [ "python", "./main.py" ]
