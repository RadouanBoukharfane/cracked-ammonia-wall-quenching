# Phase B0 : flamme libre instable jusqu'au régime saturé (pilote au point de Chabot)

## But (01/11/2026)

Tester tôt le cœur du pari, sans attendre la carte complète de la Phase A : au point de Chabot (γ = 0,43,
φ = 0,5, 10 atm, 750 K), la flamme instable saturée a-t-elle une réactivité accrue (I₀ = (s_c/s_L)/(A/L_y) > 1,
Schneider) et une variation de mélange local (φ_loc, γ_loc) importante ? Ce run fournit ensuite les états
initiaux du HOQ 2D instable (Phase B), à comparer au HOQ 1D, polluants compris.

## Configuration (`make_b0.py`)

- x périodique, L_x = 100 l_F = 11,87 mm ; y entrée → sortie, L_y = 100 l_F.
- Base 512 × 512 (l_F/5,12), 2 niveaux d'AMR dynamique sur ∇T (saut > 0,02 (T_b − T_u) = 22 K entre mailles
  voisines), maille fine 5,8 µm = l_F/20,5 (Schneider : δ_T/20).
- Perturbation initiale : 50 modes (λ de 2 à 100 l_F), 10⁻² l_F par mode, phases tirées avec la graine SEED + 1.
- Contrôle actif PeleLMeX (`active_control.*`) : bord d'attaque de l'isotherme T_track = 1195 K maintenu à
  y = 30 l_F ; τ = τ_F ; changeMax = 0,02 s_L ; pseudo-gravité. Le fichier problème (`template/`) fait suivre
  la vitesse d'entrée au contrôleur (`FCData.ctrl_V_in + (t − ctrl_tBase) ctrl_dV`), comme triple-points.
- MA + Soret, mécanisme Mei 31 espèces (10 atm), réglages de robustesse de la Phase A (δT, divU).
- 40 τ_F (15,2 ms). Plotfiles complets (espèces, taux de réaction, dégagement de chaleur) tous les τ_F ;
  checkpoints tous les 2 τ_F (états initiaux des HOQ).
- Coût estimé : ~1,7 M cellules, ~200–250 pas/τ_F → ~10 000 pas ; ~5 s/pas sur 4 nœuds → ~13 h
  (job de 24 h, reprise sur le dernier chk sinon).

## Jobs

| Job | Rôle |
|---|---|
| 8218066 | build (`build/10_000atm`, exécutable dans `bin/`) |
| 8218067 | run `runs/B0_g0.43_phi0.50_P10.0_T750` (afterok du build), 4 nœuds |

## À faire

1. Vérifier au démarrage : contrôle actif (vitesse d'entrée, position du front), raffinement (nombre de cellules).
2. Post-traitement : s_c (consommation d'O₂ intégrée, valable pour deux combustibles), longueur A de l'isoligne
   T_track, I₀(t), PDF jointe (φ_loc, γ_loc) sur le front ; critère de saturation : I₀ et A/L_y stationnaires.
3. HOQ 2D instable (Phase B) : paroi isotherme 750 K à la place de l'entrée, initialisé depuis 2–3 checkpoints
   saturés ; références HOQ 1D (PeleLMeX 1D) et HOQ 2D étroit (L_x = 2–4 l_F).

## Résultats (02/11/2026)

- Run complet : 40 τ_F, 7 204 pas, 5 h 52 sur 4 nœuds ; 0,46 à 0,65 M cellules (moins que prévu) ; aucun avertissement δT.
- Contrôle actif : V_in passe de s_L à ~2,06 s_L (plateau dès ~18 τ_F). Le bord d'attaque reste à 30 l_F ; le centre
  du front dérive vers l'aval (30 → 43 l_F) et l'épaisseur du front plissé croît (≈ 19 l_F en moyenne après 20 τ_F,
  29 l_F à 40 τ_F) : les plus grandes échelles (λ ~ L_x) se développent encore.
- `post_b0.py` (`post_b0.sbatch`, 28 processus, 22 s) → `runs/.../b0_timeseries.csv`. Pour t ≥ 20 τ_F :
  s_c/s_L = 1,81 ± 0,10 ; A/L_x = 1,39 ± 0,09 ; **I₀ = 1,30 ± 0,02** (stable dès ~10 τ_F).
  Loi de Schneider (H₂) via A0 : I₀ = 0,978 + 0,0115 Ze/Pe_cd = 1,41 → **+8 %** : la loi se transfère bien à ce
  mélange à deux combustibles (un seul point pour l'instant).
- Piège : les espèces ne sont pas dans les plotfiles (`peleLM.plot_speciesState` ignoré ; la bonne clé est
  `amr.plot_speciesState`, corrigée dans `make_b0.py`). Plotfiles avec espèces régénérés d'un pas depuis
  chk03262, chk04404, chk05588, chk06796 (~20, 26, 32, 38 τ_F) : `runs/.../species/`, job 8219139.

## Mélange local le long du front saturé (02/11/2026, `local_mixture.py`, figure `figs/local_mixture_B0_*.png`)

φ_loc = n_H/(2 n_O) (éléments) ; γ_loc = 2X_H2/(3X_NH3 + 2X_H2) (espèces). 4 états saturés (20–38 τ_F), ~3 500
points par isoligne ; résultats identiques d'un instant à l'autre.

| | φ_loc (T_track) | γ_loc (T_track) | φ_loc (T_pre = 861 K) | γ_loc (T_pre) |
|---|---|---|---|---|
| gaz frais | 0,500 | 0,430 | 0,500 | 0,430 |
| flamme plane 1D | 0,411 | 0,238 | 0,400 | 0,290 |
| flamme instable, moyenne | 0,429 | 0,246 | 0,415 | 0,299 |
| écart-type | 0,013 (3 %) | 0,035 (15 %) | 0,010 | 0,019 |
| crêtes convexes (κ l_F > 0,2) | 0,436 | 0,263 | 0,420 | 0,308 |
| creux concaves (κ l_F < −0,2) | 0,405 | 0,176 | 0,400 | 0,269 |

- γ_loc varie ~5 fois plus que φ_loc en relatif ; corrélation positive avec la courbure (0,6–0,7) : crêtes enrichies
  en H₂, creux enrichis en NH₃ (γ_loc −26 % par rapport au plan dans la zone de réaction).
- Sur l'isotherme de préchauffage, (φ_loc, γ_loc) tombent sur une seule courbe (pente dγ/dφ ≈ 2) : les deux sont
  pilotés par la même diffusion différentielle de H₂ → le mélange local est, ici, à un seul paramètre.
- Conséquence pour H1/C3b : l'ensemble de flammes 1D doit suivre ce chemin (φ et γ variant ensemble), pas φ à γ fixé
  comme pour H₂ pur. Les creux riches en NH₃ et plus pauvres sont les candidats pour le glissement de NH₃ et le N₂O
  près de la paroi : à vérifier dans les HOQ.
- Amplitude modeste (φ_loc ±3 %) comparée aux flammes H₂ pures de Schneider (TDI plus forte), cohérente avec
  I₀ = 1,30.

## Série 2 (02/11, 13 h 30) : γ et T_u

Motivation (Phase A, 750 K, 10 atm, φ = 0,5) : ω̄_max linéaire 0,43 / 0,67 / 0,73 / 0,66 pour γ = 0,2 / 0,43 / 0,75 / 1,0
(maximum vers 0,75) alors que I₀ prédit (loi de Schneider via A0) décroît : 1,57 / 1,41 / 1,19 / 1,09. Question :
lequel, croissance linéaire ou réactivité saturée, pilote l'effet sur la paroi et les émissions ?
Jobs (`jobs_set2.txt`) : B0_g0.20 (8219471), B0_g0.75 (8219491), B0_g1.00 (8219631) à 750 K / 10 atm / φ = 0,5 ;
B0_g0.43 à 500 K / 10 atm / φ = 0,5 (8219700, I₀ prédit 2,6). `amr.derive_plot_vars = mass_fractions` ajouté
(démarrage direct des HOQ depuis les plotfiles).

## Résultats série 2 (02/11, 20 h ; post_b0, moyennes t ≥ 20 τ_F, saturation vérifiée sur les 10 derniers τ_F)

| γ | I₀ mesuré | s_c/s_L | A/L_x | V_in/s_L | brosse (l_F) | I₀ Schneider (A0) | écart | ω̄_max linéaire |
|---|---|---|---|---|---|---|---|---|
| 0,20 | 1,18 ± 0,01 | 1,48 | 1,25 | 1,59 | 10 | 1,57 | **−25 %** | 0,43 |
| 0,43 | 1,30 ± 0,02 | 1,81 | 1,39 | 2,06 | 19 | 1,41 | −8 % | 0,67 |
| 0,75 | 1,30 ± 0,02 | 1,69 | 1,30 | 1,69 | 13 | 1,19 | +9 % | 0,73 |
| 1,00 | 1,10 ± 0,01 | 1,36 | 1,23 | 1,42 | 8,5 | 1,09 | +1 % | 0,66 |

- I₀ mesuré est **non monotone en γ** (maximum 1,30 entre 0,43 et 0,75), alors que la loi de Schneider évaluée via A0
  prédit une décroissance monotone : la loi ne capte pas la tendance en γ (−25 % à γ = 0,2). Ni ω̄_max seul
  (γ = 1 : ω̄ comparable à 0,43 mais I₀ = 1,10) ni la loi ne suffisent. À tester : A0 évalué côté gaz brûlés (Lehmann).
- Corrélateur : voir ../A0/JOURNAL.md (02/11 20 h 30) — Ze(1 − Le_eff,b,diff) R² = 0,86 ; Schneider R² = 0,06.
- HOQ instables lancés depuis ces états (t ≈ 26 et 38 τ_F) : voir PhaseB/jobs_set3_unstable.txt.
