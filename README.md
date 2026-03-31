# Maximum-Toroidal-Distance-Codes-DFR-Estimation
Python codes to computing the DFR of MTD codes in ML-KEM (Kyber)

Reproduce the Results (Table 2) in the paper: S. Liu and A. Sakzad, “On the Maximum Toroidal Distance Code for Lattice-Based Public-Key
Cryptography,” IEEE ISIT 2026, Accepted for Publication. https://arxiv.org/abs/2601.08452

- MTD_2D_Codes_(optimized 2D Minal Codes).py
  - The DFR of the optimized 2D Minal codes
  - Minal codes: T. B. Paiva, M. A. S. Jr, S. M. Hafiz, B. Yildiz, E. L. Cominetti, and H. S. Ogawa, “Tailorable codes for lattice-based KEMs with applications to compact ML-KEM instantiations,” in IACR Transactions on Cryptographic Hardware and Embedded Systems, Jun. 2025, p. 139–163. [Online]. Available: https://doi.org/10.46586/tches.v2025.i3.139-163
- GTD_4D_Codes_(based on the D4 lattice).py
  - The DFR of the proposed 4D GTD codes
- GTD_8D_Codes_(based on the E8 lattice).py
  - The DFR of the proposed 8D GTD codes
  - The **240 root vectors (or Voronoi vectors {v_i})** of 2*E8 lattice are used to compute the distribution **D_i=Pr{<x_i,v_i>}** and the correspoinding DFR per codeword
  - Follow the same method in: Charbel Saliba, Laura Luzzi, Cong Ling. A reconciliation approach to key generation based on Module-LWE. IEEE International Symposiumon Information Theory 2021, Jul 2021, Melbourne, Australia. https://ieeexplore.ieee.org/document/9517882
