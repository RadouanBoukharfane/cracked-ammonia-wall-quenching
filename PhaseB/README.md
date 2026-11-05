# Phase B : HOQ 2D instable contre HOQ plan (pilote au point de Chabot)

## Configuration (02/11/2026, `make_hoq.py`)

- Paroi en y = 0 : adhérente, isotherme T_w = T_u = 750 K, inerte (flux d'espèces nul) ; sortie en y = L_y ;
  x périodique. Avec Soret et paroi isotherme, PeleLMeX met `use_wbar = 0` et active la correction Soret de paroi.
- **Instable** : initialisé (`peleLM.initDataPlt`) depuis un état saturé du B0 (plotfile régénéré avec les
  fractions massiques, `../PhaseB0/runs/B0_.../species/`). Changement de repère par `patchFlowVariables`
  (template/pelelmex_prob.H) : on retire V_in(t₀) ≈ 2,06 s_L à la vitesse verticale → gaz frais au repos.
  Bord d'attaque du front à ~30 l_F de la paroi. Domaine et maillage du B0 (100 × 100 l_F, l_F/20,5) + bande
  y < 3 l_F toujours au niveau fin. 30 τ_F, 4 nœuds.
- **Plan (référence « stable », équivalente au HOQ 1D)** : domaine étroit L_x = 3,1 l_F (16 mailles de base),
  flamme 1D (PMF) dans le repère de la paroi (prob.v_shift = s_L), front à 30 l_F. 45 τ_F, 1 nœud.
- Sorties : plotfiles complets (espèces, taux de réaction, dégagement de chaleur) tous les 0,5 τ_F ;
  checkpoints tous les 5 τ_F.

## Jobs (build 8219149)

| Job | Run |
|---|---|
| 8219150 | `runs/HOQ_g0.43_phi0.50_P10.0_T750_planar` |
| 8219151 | `runs/HOQ_g0.43_phi0.50_P10.0_T750_t20.0` (état B0 à 20 τ_F) |
| 8219152 | `runs/HOQ_g0.43_phi0.50_P10.0_T750_t32.0` (état B0 à 32 τ_F) |

## Diagnostics prévus (méthode Schneider Part I)

- Flux pariétal Φ_w(x, t) = λ ∂T/∂y à la paroi ; Φ_q = max_t,x Φ_w ; normalisé par le HOQ plan.
- Pe_q = min_t (distance du max de HRR à la paroi)/δ_T, global et local (par colonne x).
- Émissions : NH₃, H₂, NO, N₂O moyennés dans 0 < y < 2 q_d après le quench (normalisation par X_NH₃ entrant et
  par l'énergie).
- Test de la loi Φ_q,2D/Φ_q,1D = −1,16 + 2,07 I₀ avec I₀ = 1,30 du B0 → attendu ≈ 1,53.

## Incident du 02/11 midi : plantage après l'extinction

- Le HOQ plan s'est arrêté à 28,8 τ_F (juste après le quench) : `checkDt` → « divU dtcell larger than actual dt ».
  Après l'extinction, l'écoulement ralentit, le pas convectif grimpe à 9·10⁻⁵ s (×15) et dépasse la contrainte
  divU au premier sursaut de dégagement de chaleur près de la paroi.
- Correctif : `amr.max_dt` = 1,5 × CFL Δx_fin / ((σ − 1) s_L) = 1,0·10⁻⁵ s ; checkpoints tous les 2 τ_F. Ajouté à
  `make_hoq.py` et aux 3 inputs. Les 2 HOQ instables ont été arrêtés préventivement et relancés (ils auraient
  planté au même endroit) ; logs dans `crash1/`.
- Relances : plan 8219211 (depuis chk01614, ~25 τ_F), t20 8219212 (depuis chk00379, 5 τ_F), t32 8219214 (depuis
  chk00000).

## Post-traitement (`post_hoq.py RUN NPROC` → `hoq_timeseries.csv`, `hoq_wallflux.npz`)

Référence plane (02/11, 12 h 30) : q_ref = ρ_u c_p s_L (T_ad − T_u) = 1,84 MW/m² ; **Φ_q* = 0,231** à 29,0 τ_F ;
**Pe_q = 1,6–1,8** (1,83 au flux maximal, 1,59 au plus près). Émissions près de la paroi (y < 2 l_F) très
dépendantes du temps après le quench : NH₃ 1,9 % → 0,04 %, N₂O 1 160 → 440 ppm, NO 1 250 → 2 000 ppm entre 28,5 et
44 τ_F. Comparaison plan/instable à faire à délais fixes après le quench (t_q + 5 et + 10 τ_F).

## Série 2 (02/11, 13 h 30)

- Références planes : g0.20 (8219735), g0.75 (8219740), g1.00 (8219743) à 750 K, g0.43 à 500 K (8219746).
- Réalisations supplémentaires au point de Chabot : t26 (8219749), t38 (8219752) → 4 réalisations.
- Nœuds lents exclus dans `make_hoq.py` (bloc h23a5-u13/u14 : ×3 sur le temps par pas).

## Incident du 02/11, 14 h 30 : HOQ « instables » partis de la flamme plane

- `peleLM.initDataPlt` est ignoré sans message : la clé est **`amr.initDataPlt`** (alors que
  `peleLM.initDataPlt_patch_flow_variables` et `peleLM.initDataPlt_reset_time` gardent le préfixe peleLM).
  Les 4 runs t20/t26/t32/t38 sont donc partis de la PMF (front plan à 59 l_F, avance à s_L) : sorties archivées
  dans `runs/HOQ_*_t*/wrong_init/` (ce sont des HOQ plans en grand domaine, non utilisés).
- Correction dans `make_hoq.py` et les inputs ; relances 8221120–8221123 (t20, t26, t32, t38). Vérifié sur
  plt00000 de t20 : lecture du plotfile du B0 sur les 3 niveaux, front plissé (y/l_F de 29,9 à 42,8, moyenne 35,7),
  gaz frais au repos en moyenne.
- Post-traitement des références planes (`post_hoq.sbatch`, job 8221118) :
  γ = 0,2 : Φ_q* = 0,290, Pe_q = 1,6–1,7 ; γ = 0,43 : 0,231, 1,6–1,8 ; γ = 0,75 : 0,209, Pe au flux maximal 1,2 ;
  γ = 1,0 : 0,210, Pe au flux maximal 0,51 (Pe_q « min » de 0,17 et 0,02 à γ = 0,75 et 1 : HRR de paroi après
  l'extinction, à redéfinir).
- Vérification du correctif (02/11, 15 h) : les 4 relances lisent le plotfile du B0 sur les 3 niveaux ; t₀ = 0
  (reset_time appliqué) ; gaz frais quasi au repos (v_y moyen 0,005–0,02 m/s, < 8 % s_L) ; le bord d'attaque passe
  de 30,2 à 29,0 l_F en 0,5 τ_F, soit ≈ 2,4 s_L, et le front moyen ≈ 1,8 s_L (cohérent avec s_c/s_L = 1,81 du B0).
  Aucune autre clé ignorée dans les logs terminés (seul l'ancien B0 listait `peleLM.plot_speciesState`, corrigé).
  Coût ≈ 280 pas/τ_F, ~4 s/pas → contact avec la paroi vers ~16 τ_F (≈ 20 h), fin des 30 τ_F en ~9–10 h.

## HOQ plan à 500 K (γ = 0,43, φ = 0,5, 10 atm ; job 8219746, post 8221236 ; 02/11, 14 h 50)

Φ_q* = 0,221 (750 K : 0,231) à 27,5 τ_F ; Pe au flux maximal 2,6 (750 K : 1,8) ; Pe_min non significatif
(HRR parasite, même problème de définition). Émissions près de la paroi (y < 2 l_F), X_NH3 entrant = 6,65 % :

| | NH₃ t_q / +5 / +10 τ_F | NO (ppm) | N₂O (ppm) | NO gaz brûlés 1D |
|---|---|---|---|---|
| 750 K | 1,48 / 0,28 / 0,09 % | 1093 / 1101 / 1544 | 1112 / 771 / 566 | 6055 ppm |
| 500 K | 2,13 / 0,62 / 0,31 % | 122 / 81 / 96 | 353 / 225 / 171 | 2026 ppm |

À 500 K : glissement de NH₃ ×1,5 à ×3,4, NO près de la paroi ÷10 (÷3 seulement en gaz brûlés 1D), N₂O devient
l'oxyde d'azote dominant près de la paroi (N₂O/NO ≈ 2–3 contre ≈ 0,4–1 à 750 K). Attention : T_u et T_w varient
ensemble (T_w = T_u) ; l'effet propre de la paroi froide n'est pas séparé.

## Série T_w (02/11, 15 h) : effet propre de la paroi froide sur NH₃ → N₂O

Question : la forte baisse de NO et la dominance de N₂O près de la paroi à 500 K viennent-elles de la paroi froide
ou des gaz frais plus froids ? HOQ plan au point de Chabot (γ = 0,43, φ = 0,5, 10 atm, T_u = 750 K) avec
T_w = 400, 500, 600, 650 K (+ référence T_w = 750 K). `make_hoq.py planar --tw TW` (nouvelle option, suffixe `_Tw`).
Jobs 8221238 (400), 8221239 (500), 8221240 (600), 8221241 (650), 1 nœud, 45 τ_F.
Vérifié à 1 τ_F : couche limite thermique en place (T de la 1re maille = T_w + 4–14 K, 748 K à 2 l_F), front
à 28,5 l_F, T_max = 1838 K, pas d'erreur.

### Résultats série T_w (post 8221520, 02/11 16 h 20)

Φ_q pris pour t > 15 τ_F (le maximum brut est à t = 0 : saut 750 K → T_w sur une demi-maille au démarrage).
« fond » = flux de conduction gaz frais → paroi 5 τ_F avant t_q. Moyennes y < 2 l_F.

| T_w (K) | t_q (τ_F) | Φ_q* (fond) | Pe au flux max | NH₃ t_q / +5 / +15 | NO t_q / +15 (ppm) | N₂O t_q / +15 (ppm) | N₂O/NO à +15 |
|---|---|---|---|---|---|---|---|
| 750 | 29,0 | 0,231 (0) | 1,83 | 1,48 / 0,28 / 0,04 % | 1093 / 1991 | 1112 / 439 | 0,22 |
| 650 | 30,0 | 0,187 (0,006) | 2,27 | 2,12 / 0,70 / 0,18 % | 675 / 1038 | 954 / 545 | 0,53 |
| 600 | 31,0 | 0,172 (0,010) | 2,56 | 2,14 / 0,87 / 0,31 % | 544 / 778 | 878 / 585 | 0,75 |
| 500 | 32,5 | 0,151 (0,026) | 3,05 | 2,43 / 1,28 / 0,65 % | 423 / 514 | 771 / 625 | 1,22 |
| 400 | 34,0 | 0,138 (0,061) | 3,44 | 2,64 / 1,64 / 1,06 % | 364 / 401 | 697 / 632 | 1,58 |

- Tendances monotones et nettes : paroi plus froide → extinction plus loin (Pe 1,8 → 3,4), flux réduit, NH₃
  résiduel ×26 à t_q + 15, NO divisé par 5, N₂O presque figé (le NH₃ et le N₂O ne sont plus consommés).
- Comparaison avec T_u = T_w = 500 K : à T_w = 500 K, T_u = 750 K donne NH₃ 0,80 % à +10 et 0,65 % à +15, N₂O/NO 1,2 ;
  T_u = 500 K donne NH₃ 0,31 % à +10, NO 96–128 ppm, N₂O/NO ≈ 1,1–1,8. La paroi froide explique
  l'essentiel du glissement de NH₃ et de la dominance de N₂O ; le niveau bas de NO à T_u = 500 K vient surtout de T_u.
- Limite : la couche limite thermique froide croît pendant les ~30 τ_F avant l'arrivée de la flamme (épaisseur
  ~ √(α t)) ; le résultat dépend donc de la distance initiale (30 l_F). À tester (distance 15 et 60 l_F).

### Série T_w, ensemble complet (02/11, 16 h 35 ; `jobs_Tw_set2.txt`)

- Test de distance (γ = 0,43) : front initial à 15 l_F (32 τ_F) et 60 l_F (78 τ_F), T_w = 400/500/600/650/750 K
  (`make_hoq.py --dist`, suffixe `_d15`, `_d60`) → l'épaisseur de la couche limite froide à l'arrivée de la flamme
  varie d'un facteur ~2 (√(t)).
- γ = 0,20, 0,75, 1,0 : T_w = 400/500/600/650 K (référence T_w = 750 K = HOQ plans série 2).
- 22 jobs 8221533–8221556 ; démarrage vérifié (front initial, T de la 1re maille = T_w + 1–5 K, pas d'erreur).
  g0.75_Tw400 tombé sur un nœud lent (slurm-compute-h21a5-u14-svn3, Advance ×6–8) : annulé, log dans
  `slow_node/`, relancé 8221598 ; nœud ajouté aux exclusions de tous les générateurs.
- `post_hoq.py` : Φ_q cherché pour t > dist/2 τ_F (le pic de démarrage à T_w < T_u est exclu).

### Résultats ensemble T_w (02/11, 17 h 40 ; `summary_tw.py` → `tw_summary.csv`, `figs/tw_series.png`)

Valeurs à t_q + 10 τ_F, y < 2 l_F. γ = 0,2 / T_w = 400 K trop court (t_q ≈ 36 τ_F) : prolongé à 55 τ_F (job 8221687).

- **Contrôle** : à T_w = T_u = 750 K, la distance initiale ne change rien (Φ_q* = 0,229 / 0,231 / 0,232,
  NH₃ 0,099 / 0,092 / 0,098 % pour d = 15 / 30 / 60 l_F).
- **Distance (paroi froide)** : effet de même signe, amplitude ±25–30 % ; à T_w = 400 K, NH₃ 0,86 / 1,13 / 1,44 %,
  N₂O/NO 1,28 / 1,63 / 1,99, Φ_q* 0,166 / 0,138 / 0,114 pour d = 15 / 30 / 60 → piloté par l'épaisseur de la
  couche limite froide à l'arrivée de la flamme (∝ √t). Paramètre de regroupement à chercher (δ_th/l_F, T à y = Pe_q).
- **NH₃ résiduel / NH₃ entrant** (d = 30) : 1,4 % (γ = 0,43 et 0,75), 0,4 % (γ = 0,2) à T_w = 750 K →
  12 % (γ = 0,43), 9 % (0,2 et 0,75) à 500 K, 17 % (0,43) à 400 K. Même ordre de grandeur pour tous les γ.
- **N₂O presque indépendant de T_w** (γ = 0,43 : 570–680 ppm ; γ = 0,75 : 710–840 ; γ = 0,2 : 190–265) : la montée
  de N₂O/NO vient de l'effondrement du NO près de la paroi froide (γ = 0,43 : 1540 → 400 ppm ; γ = 0,2 : 840 → 105),
  pas d'une production accrue de N₂O. (Corrige l'interprétation « N₂O non détruit » du paragraphe précédent.)
- **N₂O dominant (N₂O/NO > 1) seulement pour γ ≤ 0,43** et T_w ≤ 500–600 K ; à γ = 0,75 (NO 1D = 9 760 ppm),
  N₂O/NO reste < 0,35. γ = 0,2 est le plus sensible (N₂O/NO 0,23 → 2,5 entre 750 et 500 K).
- Φ_q* décroît et Pe_q croît quand T_w baisse pour tous les γ (extinction plus loin de la paroi).
- γ = 1 (sans NH₃) : NO thermique seulement (2–5 ppm), omis des panneaux d'émission.

B0 γ = 1,0 (post 8221682) : s_c/s_L = 1,36, A/L_x = 1,23, **I₀ = 1,10 ± 0,01** pour t ≥ 20 τ_F ; loi de Schneider
via A0 : 1,09 (+1 %). Deuxième point de validation après γ = 0,43 (+8 %).

### Regroupement des distances et chimie du NO (02/11, 18 h ; `coldwall_analysis.py`, `collapse_tw.py`, `no_budget.py`)

- **Instant de quench** : avec T_w < T_u, la flamme cale plusieurs τ_F AVANT le maximum du flux pariétal (retard
  1,5 τ_F à T_w = T_u → 7 τ_F à 400 K, d = 30 ; 11,7 τ_F à d = 60). Nouveau repère t_Q : ∫HRR dy < 50 % de la flamme
  libre ; Pe_Q = pic de HRR à t_Q ; émissions à t_Q + 10 τ_F (tableaux `tw_collapse.csv`). Le tableau précédent
  (t_q = flux maximal) est à remplacer par celui-ci.
- **Paramètre de regroupement** : déficit d'enthalpie de la couche froide
  H = ∫₀^{y_f−2 l_F} (ρ/ρ_u)(T_u − T) dy / ((T_ad − T_u) l_F), pris quand le bord d'attaque est à 8 l_F de la paroi.
  R² d'un ajustement quadratique unique sur d = 15/30/60 (γ = 0,43) : 0,93–0,99 (θ_w seul : 0,63–0,94 ;
  épaisseur de couche : 0,32–0,85 ; T moyenne jusqu'à y_Q : 0,83–0,98). NH₃ résiduel quasi universel en γ.
- **Bilan NO (y < 2 l_F, t_Q − 2 → t_Q + 10)**, γ = 0,43, 750 → 400 K : la couche est un puits net de NO dans tous
  les cas (le NO arrive par transport). Production ÷16 (HNO + OH ÷230, NH + O₂ ÷80, HNO + O₂ et H + NO₂ ÷10),
  consommation ÷9,5 (NH₂ + NO, les deux voies ; NH₃ + OH ÷17), transport ÷2,4. HO₂ + NO → NO₂ passe de 16 à 43 % de la
  consommation mais NO₂ est recyclé (production nette faible). Pas de nouvelle voie de destruction : le NO est bas
  parce que la flamme s'éteint plus loin et que le pool radicalaire manque dans la couche froide.
- **Bilan N₂O** : à 750 K formé sur place (NH + NO), détruit (H + N₂O, décomposition) et exporté ; à 400 K chimie
  locale négligeable, 87 % de l'inventaire apporté par transport depuis la zone de quench puis figé.
- Sous-section du papier : `../paper/sec_cold_wall.tex` (+ `figs/coldwall_*.pdf`, `coldwall_standalone.pdf`) ;
  figures générées par `paper_figs_coldwall.py`.

## HOQ instables γ = 0,43 : 4 réalisations (02/11, 20 h)

| réalisation | Φ_q* | Φ_q,2D/Φ_q,plan | t au max (τ_F) | Pe au max | Φ moyen en x max | NH₃ / NO / N₂O à t = 30 τ_F |
|---|---|---|---|---|---|---|
| t20 | 0,405 | 1,75 | 16,0 | 1,00 | 0,181 | 0,101 % / 2147 / 626 ppm |
| t26 | 0,428 | 1,85 | 24,5 | 0,95 | 0,164 | 0,094 % / 2308 / 632 ppm |
| t32 | 0,453 | 1,96 | 25,5 | 0,85 | 0,154 | 0,111 % / 2181 / 640 ppm |
| t38 | 0,413 | 1,79 | 16,5 | 1,00 | 0,149 | 0,127 % / 2010 / 671 ppm |

- Ratio moyen **1,84 ± 0,09** contre 1,53 prédit (loi Φ_q,2D/Φ_q,1D = −1,16 + 2,07 I₀, I₀ = 1,30) : +20 %.
  Pe au flux max 0,85–1,0 contre 1,83 en plan. Émissions en fin de run reproductibles (±7 %). Comparaison des émissions
  avec le plan à faire avec un t_Q fondé sur la HRR (quench étalé sur ~10 τ_F en 2D).
- Série 3 (jobs_set3_unstable.txt) : HOQ instables γ = 0,20 / 0,75 / 1,0, 2 réalisations (t26, t38), démarrage vérifié.

## HOQ instables contre paroi froide (02/11, 20 h 40 ; `jobs_unstable_coldwall.txt`)

γ = 0,43, T_u = 750 K, T_w = 500 K, depuis les mêmes états B0 que t26 et t38 (8222342, 8222343) : croise les deux
résultats (instabilités × paroi froide). Paroi froide imposée à t = 0, couche froide formée pendant l'approche
(~16 τ_F), comme dans la série plane d = 30 l_F (référence plane : HOQ_..._planar_Tw500).

## HOQ instables γ = 1,0 (02/11, 22 h 50 ; 2 réalisations t26, t38)

Φ_q* = 0,330 / 0,331 contre 0,210 en plan → **ratio 1,57** (identique entre réalisations) ; loi de Schneider
(−1,16 + 2,07 I₀, I₀ = 1,10) : 1,12 → **+40 %** (γ = 0,43 : 1,84 contre 1,54, +20 %). Le gain de flux dû aux
instabilités ne s'annule pas quand I₀ → 1 : la relation Φ_2D/Φ_1D(I₀) de Schneider sous-estime, d'autant plus que
I₀ est faible. Piste : comparer au HOQ 1D au mélange local des crêtes (argument de Schneider Part I).
Post-traitements désormais lancés en local (pas de sbatch).

## HOQ instables γ = 0,75 (02/11, 23 h ; post local) et HOQ plans « au mélange des crêtes » (`make_crest_hoq.py`)

γ = 0,75 : Φ_q* = 0,390 / 0,383 contre 0,209 en plan → **ratio 1,85** (Schneider via I₀ = 1,30 : 1,53).
Bilan Φ_2D/Φ_plan : γ = 0,43 → 1,84 ; 0,75 → 1,85 ; 1,0 → 1,57 (γ = 0,2 en cours).

Test de l'argument de Schneider (Part I) : la crête qui s'éteint se comporte comme une flamme 1D à son mélange local.
Mélange local des crêtes sur T_track (`PhaseB0/local_mixture.py`, adapté aux plotfiles de la série 2), moyenne des
crêtes convexes (κ l_F > 0,2, « cmean ») et 95ᵉ centile (« p95 ») ; inversion Newton (Cantera, 2–3 itérations) vers le
mélange frais équivalent (φ_u, γ_u) ; pseudo-B0 dans `crest/`, HOQ plan standard (`runs/HOQ_crest_*`, jobs_crest.txt).

| γ | cible cmean (φ_loc, γ_loc) | φ_u, γ_u équivalents | s_L/s_L,ref | cible p95 | φ_u, γ_u | s_L/s_L,ref |
|---|---|---|---|---|---|---|
| 0,20 | 0,465 ; 0,100 | 0,528 ; 0,217 | 1,30 | 0,468 ; 0,105 | 0,532 ; 0,223 | 1,37 |
| 0,43 | 0,436 ; 0,263 | 0,541 ; 0,453 | 1,38 | 0,443 ; 0,277 | 0,553 ; 0,467 | 1,55 |
| 0,75 | 0,401 ; 0,609 | 0,558 ; 0,768 | 1,35 | 0,411 ; 0,625 | 0,574 ; 0,777 | 1,49 |
| 1,00 | 0,381 | 0,549 | 1,22 | 0,392 | 0,564 | 1,30 |

Attendu (pente 1D de Schneider −1,16 + 2,07 s_L/s_L,ref, indicatif) : γ = 0,43 → 1,69–2,05 (mesuré 1,84) ;
0,75 → 1,63–1,92 (1,85) ; 1,0 → 1,37–1,53 (1,57) → le flux serait fixé par les crêtes les plus réactives, pas par I₀
moyen. À confirmer par les HOQ crête (flux absolu en W/m² comparé à celui du HOQ instable).

## Incident du 02–03/11 : panne d'écriture Lustre (≈ 23 h 20 → 12 h 37)

Plus aucune écriture n'a abouti pendant ~13 h (`VisMF::Write failed` dans le B0 500 K, qui a planté ; les autres jobs
sont restés bloqués en I/O jusqu'à leur limite de temps). Quota non en cause (30,8 / 54,6 To ; OST à 41 %).
Checkpoints vérifiés avec `check_chk.py` : 2 incomplets mis de côté (`bad_chk/` : crest g0.75 cmean chk00205,
A_g0.60_phi0.70_P01.0_T500 chk00447). Relances depuis le dernier chk complet (logs dans `crash_io/`) :
`jobs_restart_io.txt` (HOQ paroi froide ×2, HOQ crête ×6) ; vague 2 : 17 relances (`resubmit_unfinished.py --nice 1000
--jobs wave2_jobs.txt`, nouvelles options) + 3 jobs tués avant d'avoir écrit leur log, relancés de zéro.
B0 500 K : log archivé dans `PhaseB0/runs/B0_g0.43_phi0.50_P10.0_T500/crash_io/`, suite 8222549 (depuis chk08756).
Les jobs restés en vie (B0 tests, crête γ = 1) ont repris leurs écritures à 12 h 37.

## Rapport des flux Φ_2D/Φ_plan, 4 γ (03/11, 13 h 20)

| γ | I₀ (B0) | Φ_2D/Φ_plan | n réal. | Schneider (−1,16 + 2,07 I₀) | écart |
|---|---|---|---|---|---|
| 0,20 | 1,18 | 1,42 (1,46 ; 1,39) | 2 | 1,29 | +10 % |
| 0,43 | 1,30 | 1,84 ± 0,08 | 4 | 1,54 | +20 % |
| 0,75 | 1,30 | 1,85 (1,87 ; 1,83) | 2 | 1,53 | +21 % |
| 1,00 | 1,10 | 1,57 (1,57 ; 1,58) | 2 | 1,12 | +40 % |

Le flux suit I₀ en tendance (γ = 0,43 et 0,75, même I₀, même ratio) mais la loi sous-estime partout, de plus en plus
quand γ augmente (H₂ plus diffusif → crêtes plus enrichies). Test en cours : HOQ plans au mélange des crêtes.
Libération de nœuds (accord utilisateur) : 9 cas de la vague 2 sans chk au-delà de chk00000 (~20 min perdues chacun)
annulés et remis en file (8223729–8223737, nice 1000) ; calculs prioritaires tous démarrés.

## HOQ instables contre paroi froide (T_w = 500 K, γ = 0,43 ; 03/11, 15 h, post local)

| | Φ_q* | / plan T_w = 500 K (0,151) | / plan T_w = 750 K (0,231) | fin (t = 30 τ_F) : NH₃ / NO / N₂O |
|---|---|---|---|---|
| t26 | 0,259 | 1,72 | 1,12 | 0,87 % / 647 / 787 ppm |
| t38 | 0,238 | 1,58 | 1,03 | 1,01 % / 623 / 790 ppm |
| instable T_w = 750 K (t26 / t38) | | | | 0,09–0,13 % / 2010–2308 / 632–671 ppm |
| plan T_w = 500 K (t_Q + 10) | | | | 1,28 % / 434 / 730 ppm |

- Le gain de flux dû aux instabilités subsiste contre une paroi froide (×1,58–1,72 contre ×1,84 à T_w = T_u) ; en
  absolu, il compense la paroi froide (flux ≈ celui du plan à 750 K).
- Émissions (provisoire, délais différents) : par rapport au plan à 500 K, les instabilités réduisent le NH₃ résiduel
  (~1 % contre 1,3 %) et relèvent le NO (~635 contre 434 ppm) ; N₂O inchangé (~790) → N₂O/NO ≈ 1,25 : **N₂O reste
  dominant**. Par rapport aux flammes instables à 750 K : NH₃ ×8, NO ÷3,4, N₂O +20 %.
- ⚠ post_hoq.py : pour les HOQ instables à T_w < T_u, le pic de démarrage (t = 0) doit être exclu (filtrer t > 5 τ_F).

## Émissions à t_Q commun (03/11, 16 h ; `emissions_tq.py` → `emissions_tq.csv`, figure `paper/figs/fig_emissions.pdf`)

t_Q = HRR intégrée sur tout le domaine < 50 % de sa valeur avant contact (même définition plan / instable).
Rapports instable / plan à t_Q + 10 τ_F (y < 2 l_F) : γ = 0,2 : NH₃ ×4,9, NO ×0,8, N₂O ×1,5 (N₂O/NO 0,28 → 0,55) ;
γ = 0,43 : NH₃ ×0,9, NO ×1,5, N₂O ×1,05 ; γ = 0,75 : ×1,0 / ×1,14 / ×1,0 ; T_w = 500 K, γ = 0,43 : ×0,8 / ×1,4 / ×1,1,
N₂O/NO 1,67 → 1,26 (N₂O reste dominant). Interprétation : creux riches en NH₃ (plus marqués à γ = 0,2) → slip et N₂O ;
crêtes (Lehmann 2025) → NO. Rédigé dans `paper/main.tex` (§ Near-wall emissions, Fig. 9).

## Réalisations supplémentaires γ = 0,20 et 0,75 (03/11, 16 h 30 ; `jobs_extra_realisations.txt`)

t20 et t32 pour chaque γ (8225026–8225029) → 4 réalisations par γ. Post-traitement automatique en local :
`post_extra_realisations.sh` (attend la fin des jobs, puis post_hoq → emissions_tq → crest_at_quench ; log
`post_extra_realisations.log`). Ensuite : `cd ../paper && python make_figs.py flux emissions`, et mettre à jour
les tableaux / chiffres du § Near-wall emissions et de la Table 2 (4 réalisations pour tous les γ sauf 1,0).
