# Use the Gradle wrapper by default; override with e.g. `make GRADLE=gradle ...`
GRADLE ?= ./gradlew

# Build the plugin zip into build/distributions/
assemble:
	$(GRADLE) assemble

clean:
	rm -rf .nextflow* work build
	$(GRADLE) clean

# Unit tests (Spock)
test:
	$(GRADLE) test

# Unpack the plugin into $NXF_PLUGINS_DIR (default ~/.nextflow/plugins)
install:
	$(GRADLE) install
