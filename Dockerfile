FROM python:3.7.3-alpine3.10
COPY src /app
ENTRYPOINT ["python", "/app/app.py"]