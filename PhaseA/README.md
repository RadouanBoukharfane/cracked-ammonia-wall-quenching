# Phase A : relations de dispersion (méthode Lehmann 2024, PeleLMeX)

## Choix des cas (01/11/2026, `select_cases.py`, à partir de A0 Mei 2021 + Soret)

- Candidats : T_u ∈ {500 ; 750} K × φ ∈ {0,4 ; 0,5 ; 0,6 ; 0,7 ; 0,8} × P ∈ {1 ; 5 ; 10} atm × 8 γ = 240.
- Exclus : 28 (A0 non convergé ou S_L < 3 cm/s, donc limite d'inflammabilité), dont 25 à 500 K et φ ≤ 0,6.
  Liste : `cases/phaseA_excluded.csv`.
- **Retenus : 212 cas** (`cases/phaseA_cases.csv`, avec l_F, domaine, Δx_fin, isotherme de suivi, durée, coût).
  - Vague 1 (78 cas) : φ ∈ {0,4 ; 0,5 ; 0,6} × P ∈ {1 ; 10} atm, les deux T_u. Sert à choisir B0/B.
  - Vague 2 (134 cas) : φ ∈ {0,7 ; 0,8} et P = 5 atm.
- Validation (`cases/phaseA_validation.csv`) : Lehmann 2024 Fig. 3a / D'Alessio, X_H2,F = 0,5 sans N₂,
  φ = 0,5, T_u = 500 K, 1 et 10 bar. ⚠ Lehmann utilise le mécanisme de Zhang (30 esp.) : écart attendu avec Mei.

**Maximum d'instabilité prédit par A0** (γ minimisant Ze·(Le_eff − 1)) : γ = 0,2 à 1 atm, 0,31 à 5 atm,
0,31–0,43 à 10 atm. À 500 K et φ ≤ 0,5, le minimum tombe sur le plus petit γ inflammable (5 lignes) : ce
n'est pas un vrai minimum, la limite d'inflammabilité le tronque. À 750 K, tous les minima sont intérieurs.

## Réglages repris de Lehmann (§3.2)

- Domaine 80 l_F (x, périodique) × 20 l_F (y, entrée/sortie), base Δx = l_F/12,8 (1024 × 256),
  2 niveaux statiques pour y ∈ [7 ; 13] l_F, Δx_fin = l_F/51,2 ; front initial à y₀ = 10 l_F.
- Perturbation Al Kassar : 40 modes, A₀ = 10⁻⁶ l_F, λ ∈ [2 ; 80] l_F, phases aléatoires avec une graine fixe.
- Front suivi par l'isotherme T = T_u + 0,4 (T_b − T_u) au lieu de 1000 K (à 750 K, 1000 K serait trop près des
  gaz frais). Lehmann montre que le choix de l'isovaleur ne change pas ω.
- Durée de 8 τ_F (plateau de ω atteint en ~6 τ_F chez Lehmann).
- l_F = (T_b − T_u)/max|∇T|, avec T_b à l'équilibre.

## Coût mesuré (test de passage à l'échelle, 01/11/2026, cas Chabot γ = 0,43, φ = 0,5, 10 atm, 750 K)

1,97 M cellules (3 niveaux). Temps par pas (sans Soret) : 14,6 s (1 nœud), 8,5 s (2 nœuds), 5,6 s (4 nœuds) ;
efficacité 86 % de 1 à 2 nœuds, 76 % de 2 à 4. Avec Soret : ~7,5 s/pas sur 2 nœuds au démarrage.
≈ 415 µs·cœur par cellule et par pas (chimie à 31 espèces) : 3 à 4 fois l'estimation initiale.
Pas nécessaires ≈ 8 τ_F × 51,2 σ / CFL ≈ 2 000 (σ = 2,4) à 3 700 (σ ≈ 4,5) → 5 à 9 h sur 2 nœuds,
~500 à 900 cœur·h par cas, ~150 000 cœur·h pour les 212 cas. QoS premium-cpu : 64 nœuds max par utilisateur
(tous projets confondus), soit 32 cas simultanés à 2 nœuds ; vague 1 ≈ 1 jour de temps réel.
`make_case.py` fixe la durée du job à partir de σ (marge 30 %) ; le script reprend sur le dernier chk.

## Mise en place (01/11/2026)

- Mécanisme `Mei2021NH3` dans PelePhysics (copie de `../mechanisms/mei2021/mei2021_pele.yaml`) :
  31 espèces, 226 réactions. Écart avec le mécanisme complet < 0,1 % sur S_L, δ_T, NO, N₂O (4 flammes pauvres).
  Contient des réactions PLOG, donc CEPTR génère un mécanisme par pression : `1_000atm`, `5_000atm`, `10_000atm`,
  `0_987atm` (1 bar), `9_869atm` (10 bar). La pression est passée sur stdin (`echo P | python -m ceptr -f ...`).
- Exécutables : `bin/PeleLMeX2d.gnu.MPI.Mei2021NH3_<PDIR>.ex` (`build.sbatch`, foss/2024a ; SUNDIALS dans le
  ThirdParty partagé de PelePhysics, déjà compilé pour triple-points).
- `template/` : FlameSheet modifié (perturbation Al Kassar à `prob.nmodes` modes, phases dans `prob.phases`).
- `make_case.py NAME` → `runs/NAME/` (pmf.dat, input.2d, case.json, run.sbatch) ; `extract_front.py RUN` →
  front.npz + dispersion.csv (ω̄(k̄), DL), puis suppression des plotfiles intermédiaires.
- Pièges :
  1. PeleLMeX exige dx = dy à 10 ulp près : écrire L_y = L_x·NY/NX en pleine précision ;
  2. Soret = `transport.use_soret = 1` (et non `peleLM.use_soret`, ignoré sans message d'erreur, seulement
     dans « Unused ParmParse Variables ») ;
  3. le checkout PelePhysics partagé contient des modifications de triple-points (transport Sutherland seulement,
     `Simple` intact).

## Runs en cours

| Run | Job | Rôle |
|---|---|---|
| `runs/A_g0.43_phi0.50_P10.0_T750` | 8216968 | pilote (cas Chabot), 8 τ_F |
| `runs/val_L24_XH0.5_phi0.5_T500_P1` | 8216969 | validation Lehmann Fig. 3a, 1 bar |
| `runs/val_L24_XH0.5_phi0.5_T500_P10` | 8216970 | validation Lehmann Fig. 3a, 10 bar |

## Résultats pilote et validation (01/11/2026)

- Pilote (Chabot) : 6 h 02, 2 072 pas, 8 τ_F. Front suivi sans trou, croissance exponentielle (r² ≈ 1) pour
  k̄ < 1,7, toujours linéaire à 8 τ_F (amplitude ~10⁻³ l_F) → fenêtre de fit 1–8 τ_F.
  ω̄_max = 0,68 à k̄ ≈ 1,0 (DL : 0,52) ; ω̄/ω̄_DL ≈ 1,8 pour k̄ = 0,3–0,55 ; coupure k̄ ≈ 1,8.
- Validation (figure : `figs/validation_lehmann_fig3a.png` ; points de Lehmann numérisés :
  `cases/lehmann2024_fig3a_digitised.csv` ; tableau 1 du supplément : `cases/lehmann2024_table1.csv`).
  Chacun dans ses propres échelles : ω̄_max 0,82 contre 0,69 (1 bar) et 1,15 contre 0,96 (10 bar), soit +20 % aux
  deux pressions ; k̄_c 1,79 contre 1,57 et 1,94 contre 1,77 (+10 à 14 %). Effet de la pression reproduit
  (×1,40 contre ×1,39). Le changement d'échelle vers ses l_F, s_L ne suffit pas : bon à 1 bar (0,72), mais
  surcorrige à 10 bar (0,75 contre 0,96).
  Paramètres globaux (Mei, Cantera) : à 1 bar, Ze 9,15 contre 9,43 et Ze(Le_eff − 1) −4,46 contre −4,62 (presque
  identiques) ; à 10 bar, Ze 19,5 contre 14,9 (réponse en pression différente entre Mei et Zhang).
  **L'écart de +20 % à 1 bar n'est pas expliqué par Ze, Le_eff ou σ** ; test décisif : refaire le cas à 1 bar avec
  le mécanisme de Zhang (réf. [39] de Lehmann).
- Lehmann évalue les Le côté **gaz brûlés** (ses Le_eff 0,51 et 0,49 = nos valeurs « b » 0,513 et 0,485 ; côté
  frais : 0,578 et 0,544). La carte A0 utilisait le côté frais : à refaire côté brûlé.
- Erreurs du supplément de Lehmann : colonne k_max = copie de X_H2,F (130 lignes) ; ω̄_max = 0,86 pour
  500 K / 1 bar / φ = 0,5 / X = 0,5, contre 0,69 sur sa Fig. 3a et sa Fig. S10a (k_c = 1,58 cohérent).
- `extract_front.py` corrigé (01/11) : fenêtre arrêtée quand max A > 5·10⁻⁵ l_F et points sous 10⁻⁸ l_F exclus.
  Au-delà, le couplage non linéaire fait remonter les modes amortis (vu à 10 bar : pas de coupure). Pilote
  inchangé (ω̄_max 0,674).

## Test Zhang à 1 bar (01/11/2026)

- Mécanisme : fichiers officiels du supplément de Zhang et al., CNF 234 (2021) 111653 (`../mechanisms/zhang2021/orig/`,
  38 esp. / 263 r.) ; réduit par `../mechanisms/reduce_for_pele.py` → 31 / 248 (s_L inchangé : 30,44 contre
  30,45 cm/s). Lehmann annonce 30 / 243 : version exacte non retrouvée. Les liens GitHub fournis : `zhang.yaml`
  (EA, streamlit) = Zhang officiel sans les 3 réactions termoléculaires de Burke–Klippenstein (effet < 1 %) ;
  Tamaoki, Stagni (31/203) et un squelette (25/175) sont d'autres mécanismes.
- Flammes 1D (Cantera) contre le tableau de Lehmann (FlameMaster), 1 bar : Zhang s_L −5 %, l_F +9 %,
  Ze(Le_eff − 1) −4,63 contre −4,62 ; Mei −4 %, +10 %, −4,46. À 10 bar : Zhang −8 %, Ze 16,2 ; Mei −19 %, Ze 19,5
  (contre 14,9). Écart de code (Cantera contre FlameMaster) de 4 à 8 % sur s_L et 8 à 12 % sur l_F.
- PelePhysics `Zhang2021NH3/0_987atm` (CEPTR à 1 bar) ; `make_case.py --mech Zhang2021NH3` ; GNUmakefile `MECH=`.
- Run `runs/val_L24_XH0.5_phi0.5_T500_P1_zhang` : build 8217362, run 8217363 (afterok). Pour qu'ils passent
  devant, `nice=1000` appliqué aux jobs en attente de la vague 1 (`scontrol top` refusé).
- Attendu : si ω̄_max ≈ 0,8 (comme Mei), l'écart à Lehmann vient de la normalisation 1D, pas de la DNS ni du
  mécanisme ; si ≈ 0,69, il vient du mécanisme.
- 10 bar (lancé le 01/11 vers 18 h) : CEPTR `Zhang2021NH3/9_869atm`, build 8217992, run 8217993
  (`runs/val_L24_XH0.5_phi0.5_T500_P10_zhang`, 12 h). Flamme 1D Zhang : s_L 5,55 cm/s, l_F 348 µm (Mei : 4,88 ;
  387). Piste d'étude : sensibilité des relations de dispersion au mécanisme sous pression (Mei et Zhang
  concordent à 1 bar sur Ze(Le_eff − 1), divergent à 10 bar : Ze 19,5 contre 16,2).

### Résultat des tests Zhang (02/11, matin) — figure `figs/validation_zhang.png`

| | Mei | Zhang | Zhang en échelles de Lehmann | Lehmann |
|---|---|---|---|---|
| ω̄_max 1 bar | 0,818 | 0,761 | 0,663 | 0,690 |
| ω̄_max 10 bar | 1,147 | 1,068 | 0,879 | 0,957 |
| k̄_c 1 bar | 1,79 | 1,69 | 1,55 | 1,57–1,58 |
| k̄_c 10 bar | 1,94 | 1,89 | 1,70 | 1,77 |

- **Méthode DNS validée** : avec le mécanisme de Lehmann et ses échelles 1D, −4 % sur ω̄_max à 1 bar, −8 % à 10 bar.
- L'écart initial de +20 % = ~7 % de mécanisme (Mei contre Zhang, identique à 1 et 10 bar) + ~12 % de
  normalisation 1D (Cantera contre FlameMaster : s_L −4 à −8 %, l_F +8 à +12 %).
- Une fois normalisées, les relations de dispersion sont robustes au mécanisme (7 %) alors que Ze diffère de 20 %
  à 10 bar : la piste « incertitude cinétique » est faible en adimensionnel ; elle reste à quantifier en dimensionnel
  (τ_F Mei 7,9 ms contre Zhang 6,3 ms à 10 bar → ω dimensionnel ~15 % plus grand avec Zhang).

## Vague 1 (lancée le 01/11/2026 vers 9 h 50)

77 runs + le pilote ; job IDs dans `cases/wave1_jobs.txt`. Durées Slurm = estimation × 1,25 (6 h 40 à 12 h 20).
64 nœuds max : 32 cas tournent, les autres attendent (QOSMaxCpuPerUserLimit). Chaque job extrait le front et
écrit `dispersion.csv` à la fin. Un job arrêté par la durée se relance avec le même `sbatch run.sbatch`
(reprise sur le dernier chk, `amr.stop_time` inchangé).

### Plantages (01/11, 14 h)

- 3 runs plantés : γ = 0,43 / φ = 0,4 / 1 atm / 500 K (pas 522) et γ = 1 / φ = 0,4 / 10 atm / 500 K (pas 225) :
  `deltaT_iters not converged` (Newton sur T dans `differentialDiffusionUpdate`, 10 itérations, tol 1e-10) ;
  γ = 0 / φ = 0,6 / 1 atm / 500 K (pas 132) : `divU dtcell larger than actual dt` (`checkDt`).
- Correctif ajouté à tous les inputs et à `make_case.py` : `peleLM.deltaT_iterMax = 40`, `deltaT_tol = 1e-9`,
  `deltaT_crashIfFailing = 0` (avertissement « deltaT_iters not converged » dans pele.log : à compter),
  `peleLM.divu_dt_factor = 0.3`. Les jobs déjà en cours gardent les anciens réglages ; relance si plantage.
- Relancés : 8217364, 8217365, 8217366 (logs du plantage dans `crash1/`). `extract_front.py` ne plante plus si
  aucune fenêtre linéaire n'est encore disponible.

### Cas γ = 0 / φ = 0,4 / 1 atm / 750 K : rejeté (01/11, 17 h)

- Croissance non exponentielle, modes courts qui montent au lieu de s'amortir, ω̄(k̄) erratique (r² < 0,95).
- Modes non perturbés 41–100 à 6·10⁻⁷ l_F dès 2 τ_F, contre 10⁻¹⁰ à 2·10⁻⁹ l_F pour les 10 autres runs ; modes 32
  et 64 (découpage en boîtes) excités. Bruit numérique du même ordre que A₀.
- Écarté : résolution (zone de réaction = 12,5 mailles fines), bande raffinée (100 % du HRR dedans), pas de temps
  (limite convective, aucun avertissement). Particularités : NH₃ pur (pas de H₂ frais), flamme lente et épaisse
  (s_L 5,8 cm/s, l_F 4,5 mm, dt 3,4·10⁻⁴ s contre 3·10⁻⁵ s à γ = 0,2).
- Tests sur 2 τ_F (`diag/`, 2 nœuds) : 8217429 = A₀ 10⁻⁴ l_F ; 8217430 = ode.rtol 10⁻⁸, atol 10⁻¹⁰.
- `extract_front.py` : indicateur `quality.json` (max des modes 41–100 à 2 τ_F, alerte > 10⁻⁸ l_F) ; fin de fenêtre
  linéaire relative à A₀ (50 A₀, identique à avant pour A₀ = 10⁻⁶ l_F). `make_case.py` : options `--a0`,
  `--ode-rtol`, `--ode-atol`.

### Bruit numérique : verdict (01/11, 23 h)

- Test (a) ode.rtol 1e-8 : bruit des modes non perturbés 2,8·10⁻¹⁰ l_F (contre 6,4·10⁻⁷) → **cause = intégration
  de la chimie (CVODE rtol 1e-6) sur de grands pas de temps** ; mais ~4× plus cher.
  Test (b) A₀ = 10⁻⁴ l_F : bruit inchangé en absolu (5,9·10⁻⁷) mais 0,6 % du signal ; coût inchangé.
- Seuls les cas τ_F > ~60 ms sont touchés : γ = 0,1 / φ = 0,5 / 500 K (τ_F 92 ms, bruit 2,8·10⁻⁶) et
  γ = 0 / φ = 0,4 / 750 K (77 ms). Les 36 autres (τ_F < 30 ms) : bruit ≤ 3·10⁻⁹. Le cas γ = 0 / φ = 0,6 / 500 K
  (66 ms) tourne en CFL 0,3 : propre (10⁻¹⁰), cohérent avec un bruit croissant avec dt.
- Refaits avec A₀ = 10⁻⁴ l_F : `runs/A_g0.10_phi0.50_P01.0_T500_A1e-4` (8218688), `runs/A_g0.00_phi0.40_P01.0_T750_A1e-4`
  (8218689). Troisième cas bruité trouvé ensuite : γ = 0,1 / φ = 0,6 / 500 K (τ_F 24 ms seulement, bruit 9·10⁻⁷),
  refait aussi (`_A1e-4`, 8218697). Le bruit n'est donc pas prévisible par τ_F seul (des cas à 24–30 ms sont
  propres) : **règle = se fier à `quality.json` de chaque run et refaire avec A₀ = 10⁻⁴ l_F tout cas signalé.**
- 3 jobs « FAILED » (A_g0.43_phi0.40, A_g0.10_phi0.60, A_g0.43_phi0.60, 500 K) : calculs complets, l'extraction
  plantait sur les plotfiles `pltNNNNN.old.<id>` créés à la relance ; `extract_front.py` corrigé.

### Bilan du 02/11 matin

- Vague 1 complète (77 + pilote). Indicateur de bruit relatif : runs propres bruit/A₀ ~10⁻⁴, r² = 1, rugosité de
  ω̄(k̄) ~0,02. Reprises A₀ = 10⁻⁴ : γ = 0,1 / φ = 0,6 / 500 K acceptable (rugosité 0,11) ; γ = 0,1 / φ = 0,5 / 500 K
  et γ = 0 / φ = 0,4 / 750 K encore bruités ; nouveau cas bruité γ = 0 / φ = 0,6 / 10 atm / 750 K.
- **CFL 0,3 supprime le bruit** (γ = 0 / φ = 0,6 / 500 K : 1,6·10⁻¹⁰ avec A₀ = 10⁻⁶) → reprises `_cfl03`
  (8219000–8219002, 18 h). `make_case.py --cfl`.
- Cas suspect mais lisse : γ = 1 / φ = 0,5 / 1 atm / 750 K, ω̄_max 0,24 < DL même aux faibles k̄ ; Ze = 3,0 seulement
  (Ze(Le_eff − 1) = −1,2, σ = 2,39) ; front stable dans la bande. À comparer à la théorie de Matalon.

## Prochaines étapes

1. Suivre la vague 1 ; relancer les jobs arrêtés avant 8 τ_F.
2. Carte ω̄_max(γ, φ, P, T_u) et comparaison avec Ze·(Le_eff − 1) de A0 → choix des cas B0/B.
3. Vague 2 (134 cas).

## Reprises CFL 0,3 (02/11, 17 h 45)

- A_g0.00_phi0.40_P01.0_T750_cfl03 (8219001) terminé : bruit 3,2·10⁻¹⁰ (contre 6,4·10⁻⁷ à CFL 0,5, ÷2000),
  **ω̄_max = 0,068 à k̄ = 0,39** (DL : 0,18) contre 0,641 dans le run bruité : la valeur bruitée était fausse d'un
  facteur 10. NH₃ pur, φ = 0,4, 1 atm : quasi stable, cohérent avec Le ≈ 1. Utiliser `_cfl03` dans la carte.
- A_g0.00_phi0.60_P10.0_T750_cfl03 : bruit 9,6·10⁻¹¹, **ω̄_max = 0,097** (DL 0,23) contre 0,473 bruité (×5). NH₃ pur quasi stable.
- A_g0.10_phi0.50_P01.0_T500_cfl03 : bruit 1,5·10⁻¹⁰, **ω̄_max = 0,528 à k̄ = 0,79** (DL 0,59) contre 2,994 bruité (×6). Les 3 reprises CFL 0,3 sont faites : utiliser `_cfl03` dans la carte.

## Vague 2 (lancée le 02/11/2026 vers 20 h 45)

134 cas (φ ∈ {0,7 ; 0,8} et P = 5 atm), générés par `gen_cases.sbatch` (8222340, 8222341 ; option `make_case`
transmise). **CFL 0,3 pour les 28 cas à s_L < 0,15 m/s** (`cases/wave2_slow.txt` ; en vague 1, 4 des 5 cas bruités
avaient s_L < 0,11 m/s), CFL 0,5 pour les 106 autres (`cases/wave2_fast.txt`). Soumis avec `--nice=1000` pour laisser
passer les calculs B/B0 ; job IDs dans `cases/wave2_jobs.txt`. Contrôle : `quality.json` de chaque run.
