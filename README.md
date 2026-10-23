# Simulations thermodynamiques — Python / Streamlit

Application pour les exercices propane / n-butane / n-pentane / n-hexane et eau / méthanol : entrées modifiables, chauffage, séparation liquide–vapeur, débits et compositions des sorties, bilans matière et énergie, diagramme binaire T-x-y.

Le moteur Python reproduit les formules PRSV2-VL de **DWSIM 10.2.11**, avec ses interactions binaires et les corrélations des fichiers originaux. Il a été comparé directement au moteur officiel sur les cas des exercices et plusieurs compositions. **DWSIM et .NET ne sont pas nécessaires pour utiliser l’application.** Les références indépendantes et les seuils numériques sont dans `VALIDATION.md`.

## Lancer sous Windows

Installer Python 3.11 ou 3.12. Dans VS Code, ouvrir ce dossier puis Terminal → New Terminal et saisir :

```powershell
.\lancer.bat
```

La première exécution installe les bibliothèques (connexion Internet requise). Ouvrir l’adresse affichée dans le terminal, généralement http://localhost:8501. Garder le terminal ouvert ; Ctrl+C arrête l’application.

Alternative manuelle :

```powershell
py -3 -m venv .venv
.venv\Scripts\python.exe -m pip install -r requirements.txt
.venv\Scripts\python.exe -m streamlit run app.py
```

## Exercices et unités

- Hydrocarbures : 100 kmol/h, fractions molaires 0,1 / 0,2 / 0,3 / 0,4 ; 100 psia = 6,894757 bar absolus ; sortie 200 °F = 93,333333 °C. L’entrée proposée à 25 °C reste modifiable.
- Eau/méthanol : 50 mol/h, fractions molaires eau 0,6 / méthanol 0,4 ; entrée 25 °C ; sortie 80 °C ; pression proposée 1,01325 bar absolu.
- Les deux cas « référence DWSIM » reprennent les entrées enregistrées. Le fichier eau/méthanol contient 50,004 mol/h et une température encore à 25 °C, distincts de l’énoncé chauffé à 80 °C.
- Débits : kmol/h, mol/h, mol/s ou kg/s. Changer d’unité convertit automatiquement la valeur. Les compositions sont molaires ; une somme différente de 1 est refusée.
- Q > 0 signifie chaleur reçue. H entrée + Q = H liquide + H vapeur. Les pressions sont absolues. La référence d’enthalpie idéale est 298,15 K.

Le chauffage impose la température de sortie. Le séparateur travaille à la même température et pression, sans perte de charge ni chaleur supplémentaire. Une phase absente a un débit nul et aucune composition affichée.

## Validation et fichiers d’origine

`data/` contient les deux fichiers XML originaux et les résultats indépendants obtenus avec le moteur officiel DWSIM. L’application calcule ses résultats ; elle ne les récupère pas dans ces références. L’onglet Références permet de vérifier la concordance.

Le fichier eau/méthanol contient des liens vers un paquet thermodynamique absent. La validation utilise le seul paquet PRSV2-VL réellement sauvegardé, sans modifier l’original. Ses anciennes enthalpies enregistrées ne constituent donc pas la référence du recalcul cohérent.

Le modèle énergétique réellement appelé par PRSV2-VL est H_PR_MIX, avec le Cp idéal des fichiers. Le libellé LeeKesler présent dans le XML est un réglage hérité qui ne décrit pas cette méthode effectivement exécutée.

Les calculs couvrent les deux mélanges des exercices en régime permanent, avec une phase liquide et une phase vapeur. LLE, solides et autres équipements ne sont pas modélisés. Les tests ne constituent pas une validation expérimentale ni une certification de tous les états thermodynamiques possibles. Les échecs de convergence sont signalés et les résultats incohérents sont refusés.

## Tests

```powershell
.venv\Scripts\python.exe -m pytest tests -q
```

## GitHub et licence

Le `.gitignore` exclut l’environnement Python, les caches et les secrets. Publier les sources, `data/`, `tests/`, `README.md`, `VALIDATION.md` et `LICENSE`.

Le code adapte des équations de DWSIM, Copyright Daniel Wagner O. de Medeiros, sous **GPL-3.0-or-later** ; conserver les notices et la licence lors de la redistribution. Voir `NOTICE.md`.
