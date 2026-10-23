# Validation contre DWSIM 10.2.11

Validation du 23 octobre 2026. Le moteur portable officiel DWSIM 10.2.11.0 a été chargé avec Python.NET pour produire les références. Le moteur de l’application est une transcription Python autonome des équations PRSV2-VL ; Python.NET et DWSIM ne sont utilisés que pour la vérification indépendante.

## Couverture

- 14 états PT : les deux mélanges, plusieurs compositions, phases liquide seule, vapeur seule et équilibre biphasique ; températures 298,15 à 370 K et pressions de 101325 à 689475,729 Pa.
- 16 comparaisons de propriétés de phase : fugacités et enthalpies aux compositions de référence, phase liquide et vapeur.
- 16 températures de bulle/rosée eau/méthanol à 101325 Pa, y compris les composés purs.
- Diagramme complet : 51 points convergés à 101325 Pa. Robustesse supplémentaire : 11 points à 0,5 / 2 / 5 / 10 bar absolus.
- Tests de conservation matière/énergie, chauffage, proportionnalité au débit, composition invalide, quatre formulaires Streamlit, conversion des unités de débit et bouton du diagramme.

## Écarts maximaux mesurés

| Grandeur | Écart absolu maximal | Seuil du test |
|---|---:|---:|
| Coefficient de fugacité | 1.066e-14 | tolérance relative 10⁻¹⁰ + absolue 10⁻¹² |
| Enthalpie d’une phase (J/mol) | 0.00649377 | 0,01 |
| Fraction vapeur molaire | 5.25845e-07 | 10⁻⁶ |
| Enthalpie du mélange après flash (J/mol) | 0.0209052 | 0,05 |
| Température bulle/rosée (K) | 6.43361e-07 | 10⁻⁵ |

Les références de flash utilisent des tolérances PT demandées de 10⁻¹². Les solveurs ne s’arrêtent pas exactement au même point ; un essai supplémentaire à 10⁻¹⁸ ne modifie pas les résultats de référence. Les écarts résiduels de fraction de phase sont inférieurs à 10⁻⁶. Les très petits écarts d’enthalpie de phase viennent de l’intégration du Cp. Les valeurs de référence sont conservées dans `data/dwsim_validation.json`.

## Exercice eau/méthanol à 80 °C

Entrée : 25 °C, 101325 Pa, 50 mol/h, z eau = 0,6 et z méthanol = 0,4. Chauffage jusqu’à 80 °C, puis séparation à même pression.

| Résultat Python | Valeur |
|---|---:|
| Fraction vapeur molaire | 0.288262448 |
| Débit vapeur (mol/h) | 14.41312238 |
| Débit liquide (mol/h) | 35.58687762 |
| x eau dans le liquide | 0.700472233 |
| y eau dans la vapeur | 0.351927933 |
| Puissance de chauffage (kW) | 0.219997150 |

Puissance obtenue à partir des enthalpies du moteur officiel : 0.219997302 kW. L’écart est 1.518e-07 kW.

## Pourquoi certaines captures anciennes diffèrent

Le fichier eau/méthanol affecte aux objets le paquet PP-b14f2c81-1ace-4c17-95a2-3488cea93e3b, qui n’existe pas parmi les paquets sauvegardés. Le seul paquet présent est PRSV2-VL, PP-cd3528fc-cef1-4c91-9a43-2fbd94ca4230. Le recalcul cohérent rattache les objets à ce paquet uniquement dans la copie de validation. Les fichiers originaux restent inchangés.

À 25 °C, l’ancienne enthalpie globale sauvegardée est −43325,284386835 J/mol ; le paquet PRSV2-VL réellement sauvegardé calcule −42762,231476707 J/mol, valeur reproduite par Python. Il serait incorrect de forcer les nouvelles équations à restituer l’ancien résultat associé à une configuration manquante.

Le fichier hydrocarbures utilise des tolérances de flash de 10⁻⁴ ; sa fraction vapeur sauvegardée (0,2196244019) diffère légèrement du calcul plus convergé (0,2196005123 en Python ; 0,2196005712 avec DWSIM plus strict).

## Modèle énergétique vérifié

La classe PRSV2VLPropertyPackage de DWSIM 10.2.11 appelle H_PR_MIX. Le libellé LeeKesler hérité du XML ne sélectionne pas Lee–Kesler pour cette méthode. La transcription conserve la formule d’enthalpie de DWSIM et intègre la corrélation ChemSep Cp n°16 des fichiers, avec la même référence à 298,15 K. Les paramètres κ sont recherchés par les noms enregistrés, y compris les replis PR du moteur officiel pour les noms N-butane / N-pentane / N-hexane. Les paramètres eau/méthanol sont directionnels, comme dans cette version de DWSIM.

Cette concordance valide les exercices et la grille indiquée ; elle ne prouve pas une identité bit à bit, une exactitude expérimentale ou la validité pour des composants et des opérations non implémentés.
