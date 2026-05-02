# ---------------------------------------------------------------------------
# weightmatrix – Docker image for running and testing the package
# ---------------------------------------------------------------------------
# Build:  docker build -t weightmatrix .
# Test:   docker run --rm weightmatrix
# Shell:  docker run --rm -it weightmatrix bash
# ---------------------------------------------------------------------------

FROM python:3.11-slim

LABEL maintainer="weightmatrix"
LABEL description="Weight matrix formalism for parametric detector track-resolution prediction"

# Set working directory
WORKDIR /app

# Install build tools and copy dependency list first for layer caching
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of the source
COPY . .

# Install the package in editable mode so tests can import it
RUN pip install --no-cache-dir -e ".[dev]"

# Default command: run the test suite
CMD ["pytest", "--tb=short", "-v"]
