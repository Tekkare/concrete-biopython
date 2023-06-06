# DNA Encryption and Analysis with TFHE / concrete-biopython library

* **Overview:** This project aims to demonstrate the use of Zama Concrete library for FHE in bioinformatics, with an emphasis on DNA sequence encryption and analysis. The project will focus on creating simple BioPython DNA functions with FHE, conducting basic DNA analysis and treatment with FHE, and discussing the potential for implementing alignment algorithms with FHE.  The project's ultimate goal is to produce tooling that could serve as a foundation for a future "concrete-biopython" TFHE-enabled BioPython library, while showcasing the potential of encrypted genomic data analysis while ensuring data privacy.
* **Related links and reference:**:
  - [BioPython Library](https://biopython.org/)
  - [Smith-Waterman](https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm)
  - [BWA-MEM Algorithm](https://github.com/lh3/bwa)

## Description:
###  DNA Encryption and Analysis with Privacy-Preserving FHE Encryption using Zama Concrete Library
### Summary:
 The primary objectives include implementing simple BioPython DNA functions with FHE, conducting basic DNA analysis and treatment, and creating alignment algorithms starting from Smith-Waterman and suggestions for BWA implementation. The final goal is to contribute towards creating a TFHE-enabled BioPython library, which offers an extensive set of tools for encrypted genomic data analysis.

## Milestones:

### 1. Simple BioPython DNA functions converted for TFHE
* **Tasks:**
  * Identify simple BioPython DNA functions such as reverse complement, transcription, and translation.
  * Adapt these functions to work with FHE using the Concrete library.
  * Package these functions in a library-like format for easy reuse.
* **Deliverables:** A library of FHE-enabled BioPython DNA functions, along with a notebook demonstrating their use and performance.

### 2. DNA Analysis & Treatment
* **Tasks:**
  * Implement basic DNA analysis functions such as Hamming distance calculation, sequence sorting, and origin of replication analysis.
  * Encrypt these functions using the Concrete library to ensure privacy-preserving analysis.
* **Deliverables:** A notebook showcasing the implementation and performance of the FHE-enabled DNA analysis functions.

### 3. Alignment Algorithms
* **Tasks:**
  * Implement the Smith-Waterman alignment algorithm for sequence alignment.
  * Discuss potential approaches for BWA algorithm and other advanced alignment algorithms with FHE. !! implementation is not in the scope as many challenges are expected.
  
* **Deliverables:** A notebook detailing the implementation of the Smith-Waterman alignment algorithm with FHE, along with a discussion on potential approaches for adapting the BWA algorithm with FHE.

### 4. Documentation and Tutorial
* **Tasks:**
  * Provide a summary of all issues encountered with Concrete for improvement and potential path forward for concrete-biopython.
  * Create a comprehensive jupyter notebook detailing the project, methodologies, techniques, and results.
  * Provide clear explanations and code samples for developers to follow as a tutorial.
* **Deliverables:** A jupyter notebook, including code snippets and visualizations, illustrating the application of FHE encryption using Concrete for DNA sequence analysis.
