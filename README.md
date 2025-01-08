Projects and Programs

This repository showcases a collection of scientific tools and protocols developed to streamline data analysis and experimental design in molecular and cellular biology research. Below is a summary of each project:

1. PLISH Probe Designer

Description: A tool for designing proximity ligation in situ hybridization (PLISH) probes to detect RNA molecules within tissues. The program supports the PLISH protocol as described in studies by the Harbury/Desai labs (PubMed, eLife, Patent).

Features:

Facilitates the design of RNA-specific probes compatible with the PLISH protocol.

Outputs probe sequences optimized for in situ hybridization experiments.

Usage:

Input target RNA sequence and parameters for probe design.

Generate and visualize probe configurations.

2. Image Analysis Protocol: Alveolar Airspace Measurements

Description: A protocol for measuring alveolar features using images from HPA lectin-stained tissue sections.

Features:

Calculates alveolar airspace dimensions and related metrics.

Tailored for applications in pulmonary biology.

Usage:

Upload HPA lectin-stained images.

Execute analysis to obtain quantitative measurements.

3. Pixel Counter

Description: A program for quantifying pixels that exceed a specified intensity threshold within an image.

Features:

Works with multiple image formats.

Outputs total and relative pixel counts for the selected threshold.

Usage:

Provide an image and a threshold value.

Analyze and retrieve pixel count data.

4. Percent of Total Counter

Description: Calculates the ratio of one feature within a larger feature in an image, such as the percentage of lineage labels within membranes of specific cell types.

Features:

User-friendly interface for defining x and y features.

Generates coverage ratios with visual overlays.

Usage:

Input images and specify the features for comparison.

Retrieve percentage calculations with graphical representations.

5. Cell Counter

Description: A tool for counting nuclei in an image to quantify the number of cells present.

Features:

Automated detection of nuclei.

Provides count metrics with visualization.

Usage:

Upload images containing nuclei.

Execute the counting algorithm and retrieve cell count data.

Installation and Dependencies

Each program is implemented in Python and requires the following dependencies:

Core Libraries:

numpy: Numerical operations, array manipulations.

pandas: Data organization, manipulation, and analysis.

os: File and directory management.

warnings: Suppressing warnings for smoother outputs.

Image Analysis Libraries:

matplotlib: Data visualization and plotting results.

opencv-python: Image processing tasks like reading, writing, and manipulating images.

scipy: Advanced image filtering, morphological operations, and feature extraction.

scikit-image: Modules for segmentation, feature detection, and measurement (e.g., regionprops, label, perimeter, clear_border, watershed).

Pillow: Basic image manipulation and annotation.

Molecular Biology Libraries:

BioPython: Sequence data retrieval and handling (e.g., Entrez, SeqIO, Seq, Blast).

melting: Calculating melting temperatures of sequences.

Utility Libraries:

collections.Counter: Counting occurrences of elements in data.

difflib.SequenceMatcher: Comparing sequences or regions for similarity.

re: Pattern matching for sequence analysis.

Install dependencies with:

pip install -r requirements.txt

License

This repository is open-source and available under the MIT License.

Contributions

Contributions are welcome! Please submit issues or pull requests to improve these tools.
