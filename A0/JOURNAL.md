# Phase A0 : journal

## 01/11/2026 : mise en place et lancement (mécanisme Mei 2021)

**Mécanisme.** Mei, Zhang, Shi, Xi, Li, CNF 231 (2021) 111472. Fichiers officiels du supplément Elsevier
(`mmc2.zip`) dans `../mechanisms/mei2021/orig/` : SMM#3 (mécanisme), SMM#4 (thermo), SMM#5 (transport), SMM#2
(S_L expérimentales). 40 espèces, 257 réactions (et non 38/265 comme écrit dans la v1 du proposal, corrigé).
Conversion `ck2yaml` → `mei2021.yaml` :
- espèce fictive FN2 retirée (inerte « thermal effect » de l'article, seulement collisionneur) → 39 espèces ;
- `explicit-third-body-duplicates: mark-duplicate` (réactions de Burke–Klippenstein H+O₂+R, comme Chemkin) ;
- transport HNOH en double (même valeur), ignoré.
HONO, NO₂, HNO₂, H₂NO présents (critère (i) du §4.2 rempli).
Pour PeleLMeX, AR, HE, OH*, CO, CO₂, HCO, CH₄, C₂H₆ pourront être retirés (≈ 31 espèces) : à faire en A.

**Environnement.** `../.venv` (Python 3.12, Cantera 3.2.0, numpy, scipy, pandas, matplotlib, openpyxl, pypdf).

**Code.** `a0lib.py` (définitions, voir docstring), `run_a0.py` (run / validate / collect), `a0.sbatch`.
- T_b = T_ad à l'équilibre (HP) : la fin du domaine n'est pas à l'équilibre pour NH₃ (2106 K contre 2150 K).
- Ze : méthode Lehmann (ΔY_N2 = ±0,003, φ fixé) = valeur principale ; méthode ΔT_u = ±10 K gardée en contrôle
  (écart notable : 9,9 contre 14,6 à γ = 0,4, φ = 1, 1 atm, 298 K).
- Le des espèces dans les gaz frais (u) et brûlés (b) ; Le_fuel molaire / diffusivités / dégagement de chaleur.
- Pe_cd = max|ṁ dY/dx| / max|d j/dx|, avec d j/dx = ω̇W − ṁ dY/dx (bilan 1D stationnaire exact).
- I0_schneider = 0,978 + 0,0115 Ze/Pe_cd : prédiction de la loi H₂, pour information.

**Grille.** γ ∈ {0 ; 0,1 ; 0,2 ; 0,31 ; 0,43 ; 0,6 ; 0,75 ; 1} × φ ∈ {0,4 ; 0,5 ; 0,6 ; 0,7 ; 0,8 ; 1,0} ×
P ∈ {1 ; 5 ; 10} atm × T_u ∈ {300 ; 500 ; 750} K = 432 cas, avec et sans Soret.

**Tests.** γ = 0,4, φ = 1, 1 atm, 298 K : S_L = 36,9 cm/s contre 35,7 ± 1,1 mesuré (+3,4 %), 22 s.
γ = 0,43, φ = 0,5, 10 atm, 750 K (Chabot) : S_L = 31,2 cm/s, δ_T = 119 µm, Ze = 10,4, Le_eff(u, diff) = 0,55,
Pe_cd,H2 = 0,28, I₀ prédit 1,41.

**Jobs.** 8216939 (A0_soret + validation) et 8216940 (A0_nosoret), 1 nœud chacun.
Sorties : `results/mei2021[_nosoret]/{cases/*.json, profiles/*.csv, a0_summary.csv}`,
`results/validation_mei2021.csv`. Relancer le même sbatch reprend là où il s'est arrêté.

**Reste à faire en A0.** Markstein (flammes étirées), dépouillement et figures, mécanismes de contrôle (Jiang,
Stagni), choix des cas de la Phase A.

## 01/11/2026 (suite) : bilan des runs

- Jobs 8216939/8216940 annulés après ~18 min : les 7 derniers cas (tous à 300 K) bloquaient le solveur jusqu'au
  timeout de 1 h. Relance 8216951 (`a0_retry.sbatch`, `run_a0.py retry`) : continuation en φ depuis le cas convergé
  voisin, 10 min max par cas. 15 cas récupérés par grille.
- Bilan : 423/432 cas OK avec Soret, 426/432 sans. Les échecs restants sont tous à φ = 0,4 (un à 0,5), γ ≤ 0,31,
  surtout à 300 K, donc à la limite d'inflammabilité : on les abandonne. 76 à 80 cas convergés ont S_L < 3 cm/s
  (marqueurs vides sur les figures), à exclure de la Phase A.
- Validation Mei (38 points, 298 K) : écart moyen 4,0 % (Soret) / 4,4 % (sans) ; 26/38 dans la barre d'erreur ;
  pire écart en riche à 10 atm (+18 % à φ = 1,4).
- **Ze·(Le_eff − 1) (Le_eff diffusivités, gaz frais) est non monotone en γ** à 500 et 750 K pour φ ≤ 0,8 ;
  minimum entre γ = 0,2 et 0,43, qui se déplace vers les γ élevés quand P augmente (750 K : 0,2 à 1 atm,
  0,31 à 5 atm, 0,31–0,43 à 10 atm). À φ = 1, Le_eff − 1 change de signe vers γ ≈ 0,5–0,7.
  Minimum plus profond à 500 K (jusqu'à −13) qu'à 750 K (jusqu'à −6,6) : confirme R1, les deux T_u restent utiles.
  À 300 K le minimum tombe souvent au bord de la zone inflammable : peu exploitable.
- Ze (variation de T_u) / Ze (dilution N₂) = 1,5 en moyenne (0,86 à 3,2) : la méthode compte, comme l'écrit Lehmann.
- Soret sur S_L : quelques %, jusqu'à −10 % à γ = 1 ; positif (≤ +3 %) vers γ ≈ 0,3 à haute P.
- Figures : `results/mei2021/figs/` (sL, delta_T, tau_F, Ze, Le_eff, ZeLe, I0, Le_eff_definitions,
  soret_effect_sL) et `results/mei2021_nosoret/figs/`.

## 02/11/2026, 20 h 30 : quel corrélateur pour I₀ (mélange à deux combustibles) ?

I₀ saturé mesuré (B0, φ = 0,5, 10 atm, 750 K) : 1,18 / 1,30 / 1,30 / 1,10 pour γ = 0,2 / 0,43 / 0,75 / 1.
Ajustement linéaire I₀ = a + b X sur ces 4 points (R²) :
- Ze/Pe_cd,H2 (convention Schneider) : 0,06 ; la loi H₂ (a, b de Schneider) donne −25 % à +9 %.
- Ze/Pe_cd mixte (pondéré X_H2,F / X_NH3,F) : 0,05.
- Ze·(1 − Le_eff), Le_eff molaire ou HRR (u ou b) : ≤ 0,05.
- Ze·(1 − Le_eff), diffusivités, gaz frais : 0,65 ; **gaz brûlés : 0,86** (résidus ≤ 4 %) → Lehmann a raison
  d'évaluer Le côté gaz brûlés. Mais 4 points pour 2 coefficients : indice seulement.
Test hors échantillon, γ = 0,43 à 500 K (B0 en cours, post à 15 τ_F, NON saturé : brosse 21 → 39 l_F, V_in 3,4–5 s_L) :
I₀ ≈ 1,8 (1,71–1,89 entre 9 et 15 τ_F). Prédictions : Schneider 2,61 (+45 %), Ze(1 − Le_eff,b,diff) linéaire 2,46,
Ze(1 − Le_eff,u,diff) 2,08, Pe_cd mixte 1,37. Aucune ne tient à l'extrapolation (X de 4–5 à 11) : la réponse de I₀
est sous-linéaire aux fortes instabilités, ou le domaine 100 l_F limite les grandes échelles à 500 K. À conclure en fin
de run (40 τ_F) ; envisager un domaine 200 l_F à 500 K.

### Suite (21 h) : extrapolation ou domaine ?

- λ_max linéaire (Phase A) ≈ 6 l_F à 500 et 750 K → L_x = 100 l_F ≈ 17 λ_max : l'échelle linéaire n'est pas limitée.
- Spectre du front (bord d'attaque, `PhaseB0/front_spectrum.py`) : en régime saturé, les modes n = 1–2 (λ = L_x, L_x/2)
  portent 75–89 % de la variance pour γ = 0,43 et 0,75 à 750 K (≈ 25 % à γ = 0,2, ≈ 40 % à γ = 1), et déjà 50–77 %
  à 500 K vers 8–14 τ_F ; rms du front en croissance continue. Les grandes cellules atteignent la taille du domaine.
- Mais Schneider Part I, Fig. 11 : I₀ quasi indépendant de la taille du domaine (1,35–1,45 pour λ = 16–100 δ_T) alors
  que s_c/s_L et A/L_y croissent → le domaine gonfle s_c et A, pas I₀. Hypothèse « domaine » peu probable pour I₀.
- Schneider Part II, Fig. 12b : la loi est calée sur Ze/Pe_cd ∈ [0, ~80], I₀ ∈ [1, ~2]. Nos points 750 K (X ≤ 52) sont
  dans la gamme ; le point 500 K (X = 142) est une extrapolation à ~1,8× la borne → hypothèse retenue : réponse
  sous-linéaire de I₀ au-delà de X ~ 80, que la loi linéaire ne capte pas.
- Tests lancés (`PhaseB0/jobs_I0_test.txt`) : B0 γ = 0,43, 500 K, L_x = 50 l_F (8222298, test domaine sur notre cas) ;
  B0 γ = 0,43, φ = 0,4, 750 K, X = 65 (8222299) ; B0 γ = 0,60, φ = 0,5, 500 K, X = 100 (8222300). 30 τ_F chacun.
  `make_b0.py --lx/--nx` (dx et L_y = 100 l_F conservés, suffixe `_Lx`).

### 03/11, 15 h : tests de I₀ (domaine, X = 65)

- **Domaine** : B0 γ = 0,43, 500 K, L_x = 50 l_F : I₀ = 1,77–1,82 (15–24 τ_F) contre 1,74–1,76 pour L_x = 100 l_F
  (s_c/s_L 3,1 contre 3,9 ; A/L_x 1,8 contre 2,2) → **I₀ indépendant du domaine** (à ~3 %), comme chez Schneider ;
  le domaine agit sur s_c et A, pas sur leur rapport. Hypothèse « domaine » écartée.
- **X = 65** (γ = 0,43, φ = 0,4, 750 K, dans la gamme de calage de Schneider) : I₀ = 1,46–1,48 (saturé dès 10 τ_F)
  contre 1,73 prédit (−15 %). La loi H₂ surestime donc aussi **dans** sa gamme : l'hypothèse « extrapolation » ne
  suffit pas. L'ajustement linéaire en Ze(1 − Le_eff,b) (fait sur les 4 γ) prédit aussi 1,73.
- **Bilan sur 6 points** (4 γ à 750 K, φ = 0,4, 500 K) : Ze(1 − Le_eff,b) = 4,19 / 4,56 / 4,83 / 5,26 / 7,23 / 10,81 →
  I₀ = 1,10 / 1,18 / 1,30 / 1,30 / 1,47 / 1,75 : **classement monotone** (Ze/Pe_cd ne l'est pas), croissance
  **concave** (pente ~0,2 par unité entre 4 et 5, ~0,08 au-delà). Ze(1 − Le_eff,b) est le bon paramètre d'ordre ; la loi
  est sous-linéaire. Point X = 100 (γ = 0,6, 500 K) en cours pour compléter.

### 03/11, 20 h 15 : test hors échantillon de la loi de réactivité

γ = 0,6, T_u = 500 K (X = Ze(1 − Le_eff,b) = 10,29 ; Ze/Pe_cd = 100) : I₀ = 1,67–1,68 (14–23 τ_F). Loi
I₀ = 1 + 0,63 ln(X/3,33) ajustée sur 6 flammes : 1,71 (+2 %) ; loi H₂ : 2,13 (+27 %). Ajouté à fig_I0 (étoile).
Vague 2, φ = 0,5, 5 atm, 750 K : ω̄_max = 0,43 / 0,57 / 0,62 / 0,65 / 0,65 / 0,51 (γ = 0,2 / 0,31 / 0,43 / 0,6 / 0,75 / 1)
contre 0,43 / 0,58 / 0,67 / 0,71 / 0,73 / 0,66 à 10 atm : même forme non monotone, plus faible à 5 atm, chute plus marquée
à γ = 1.
