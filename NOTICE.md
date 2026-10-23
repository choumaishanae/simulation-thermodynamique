# Sources et licence

Le module `dwsim_model.py` adapte les équations publiques DWSIM :

- Copyright 2012–2024 Daniel Wagner O. de Medeiros ; licence GPL-3.0-or-later.
- Version de référence : **DWSIM 10.2.11**, tag `v10.2.11` du dépôt DanWBR/dwsim10.
- [PRSV2-VL : fugacités, règles de mélange et enthalpies](https://github.com/DanWBR/dwsim10/blob/v10.2.11/engine/DWSIM.Thermodynamics/PropertyPackages/Models/PRSV2-VL.vb).
- [Paquet et orientation des interactions](https://github.com/DanWBR/dwsim10/blob/v10.2.11/engine/DWSIM.Thermodynamics/PropertyPackages/PengRobinsonStryjekVera2VL.vb).
- [Paramètres purs](https://github.com/DanWBR/dwsim10/blob/v10.2.11/engine/DWSIM.Thermodynamics/Assets/prsv2.dat).
- [Interactions binaires](https://github.com/DanWBR/dwsim10/blob/v10.2.11/engine/DWSIM.Thermodynamics/Assets/prsv2_ip_vl.dat).
- [Corrélations de propriétés](https://github.com/DanWBR/dwsim10/blob/v10.2.11/engine/DWSIM.Thermodynamics/PropertyPackages/Base/PropertyPackage.vb).

Adaptations : transcription en Python, solveur PT par Rachford–Rice et substitutions successives, résolution du diagramme binaire par moindres carrés, intégration du Cp avec SciPy, contrôles de convergence et interface Streamlit. Aucun moteur DWSIM binaire n’est redistribué dans ce dossier.

Le texte intégral de la GPL version 3 figure dans `LICENSE`. Conserver cette licence et les notices de provenance lors de la redistribution du code adapté. Les bibliothèques tierces conservent leurs licences respectives.
