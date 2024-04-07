This python code computes the numerical spectra for two-neutrino double-beta decays (2v2b) without neglecting neutrino masses, as needed for searching for massive sterile neutrinos and other exotic fermions, proposed in [arXiv:2012.09281.](https://arxiv.org/abs/2012.09281)

The main file is ``get_spectrum_2v2b.py``. There you can choose between different nuclei tuning the ``mother`` parameter, the masses of the two neutrinos with ``mnu1`` and ``mnu2``, 
and which spectra you want to compute (as defined in [arxiv:1209.5722](https://arxiv.org/abs/1209.5722)):
* Spectrum wrt summed energy: always computed, output in ``mother_sums.txt`` file.
* ``dosumonly=False`` to compute also single energy spectrum ``mother_ses.txt`` and angular correlation ``mother_cor.txt``
* ``do2D=True`` to compute 2D distributions in ``mother_2ds.txt`` (heavier and time consuming).

Run it as ``python3 get_spectrum_2v2b.py``

--- 
**Contact for additions, comments and suggestions:**
Xabier Marcano (xabier.marcano@uam.es)


--- 
**Citation info:**

Search for Light Exotic Fermions in Double-Beta Decays  
M. Agostini, E. Bossio, A. Ibarra, X. Marcano.  
[``Phys.Lett.B 815 (2021) 136127``](https://www.sciencedirect.com/science/article/pii/S0370269321000678?via%3Dihub)
[``arXiv:2012.09281``](https://arxiv.org/abs/2012.09281)

```
@article{Agostini:2020cpz,
    author = "Agostini, Matteo and Bossio, Elisabetta and Ibarra, Alejandro and Marcano, Xabier",
    title = "{Search for Light Exotic Fermions in Double-Beta Decays}",
    eprint = "2012.09281",
    archivePrefix = "arXiv",
    primaryClass = "hep-ph",
    reportNumber = "TUM-HEP 1306/20",
    doi = "10.1016/j.physletb.2021.136127",
    journal = "Phys. Lett. B",
    volume = "815",
    pages = "136127",
    year = "2021"
}
```
