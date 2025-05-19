#import "@preview/ouset:0.2.0": underset, overset

#set page(
  paper: "us-letter",
  header: align(right)[],
  numbering: "1",
  margin: (
    top: 2cm,
    bottom: 2cm,
    x: 1.9cm
  )
)
#set par(
  justify: true
)
#show link: set text(blue)

#align(center, text(18pt)[
  *Modelling Single Cell RNA-Sequencing Data*
])

== Modelling A Single Homogenous Cell Population
The current view that transcription itself is stochastic and the probe used to sample transcripts is imprecise poses a problem for inference. The implications are that:
+ Even for a homogenous cell population the true transcript counts will vary from cell to cell despite an unchanged transcriptional process.
+ Differences in sample processing can introduce systematic variation where there is none.

The first implication has been addressed by empirical observations that the distribution of transcript counts is well-explained by a negative binomial, where--- within a homogenous cell population ---for gene $g$ and cell $c$ the observed number of counts is represented by:

$ Y_(g,c) " " ~  "NegBinom"(mu_g, phi.alt_g) $

with probability mass function:

$
P_("NB")(n | mu, phi.alt) 
=
binom(n + phi.alt - 1, n) 
binom(mu , mu + phi.alt)^n 
binom(phi.alt , mu + phi.alt)^phi.alt
$

and moments:

$
bb(E)(Y) &= mu \
"Var"(Y) &= mu + frac(mu^2, phi.alt)
$

The second implication highlights that due to inefficiencies in RNA capture, sequencing, and mapping; it is unlikely that directly estimating $mu_g$ will be equivalent to estimating an absolute measure of expression. Instead, we assume the true expression parameter $eta_g$(i.e, the one we would have been able to estimate had we measured all transcripts within a cell) is scaled by a batch-specific "size factor" $zeta_b$— representing how well the sample was processed:

#align(center)[
  #block(
    fill: rgb("#C17985"),
    inset: 8pt,
    radius: 5pt,
  )[
    $
    Y_(b,g,c) " " &~  "NegBinom"(mu_(b,g), phi.alt_g) \
    mu_(b,g) &= zeta_b dot eta_g \
    $
  ]
]

With multiple batches $n_b$ > 1, one might assume this would allow us to estimate all parameters as $eta_g$ is constant across batches(fixing $g$), while $zeta_b$ is constant across genes(fixing $b$).


Unfortunately, without additional information, the model above is #text(rgb("#C17985"))[unidentifiable].

#block(
  fill: luma(235),
  inset: 8pt,
  radius: 5pt,
)[
  === A Tangent on Models and Identifiability
  - A *statistical model* $cal(M)$ is a set(or _family_) of probability distributions on a given sample space $cal(S)$. A parameterized statistical model is a set of parameterized probability distributions $cal(M) := {pi_theta : theta in Theta}$ where there exists a known mapping from the parameter space $Theta$ to the model $cal(M)$.
  - Each element in $cal(M)$ is also referred to as a *model configuration*.
  - A parameterized statistical model $cal(M)$ is *identifiable* if two equal model configurations will always have the same parameters, i.e, $pi_(theta) = pi_(theta') ==> theta = theta'$ for all $pi_theta, pi_theta' in cal(M)$. \
  
  To see this in an example, let $cal(M) := {pi_theta : theta in Theta}$ be a parameterized model where $theta := (zeta, eta, phi.alt)$ and \ $pi_theta := P_("NB")(n | zeta dot eta, phi.alt)$.
  
  Notice that if $theta = (zeta, eta, phi.alt)$ 
  and 
  $theta' = (zeta', eta', phi.alt')$ 
  are two distinct parameters such that for an arbitrary constant $k in bb(R)^+$:
  - $zeta = k zeta' $
  - $eta = frac(1,k) eta'$,
  - $phi.alt = phi.alt'$
  
  Then the two model configurations(probability distributions) $pi_theta$ and $pi_theta'$ are equal despite distinct parameters: 
  $
  pi_theta 
   &= 
  P_("NB")(n | zeta dot eta, phi.alt) \
   &=
  P_("NB")(n | (frac(k, k))zeta dot eta, phi.alt) \
   &=
  P_("NB")(n | k zeta dot frac(1, k)eta, phi.alt) \
   &=
  P_("NB")(n | zeta' dot eta', phi.alt') \
   &=
  pi_theta'
  $
  
  Implying the model is unidentifiable.
]

This tells us that estimating an absolute measure of expression is practically impossible with data from a typical single-cell experiment. Instead, we settle for the next best thing by substituting the original model with one that constrains the size factors(now denoted by the variable $s$) to sum to 1:

$
arrow(1) dot arrow(s) = sum_(b = 1)^(n_b) s_b = 1
$

Where $n_b$ denotes the number of batches.

This provides with a new model that is 1) #text(rgb("#6bb091"))[identifiable] and 2) retains relative differences in magnitude between relative expression quantities(now denoted by the variable $q$):

#align(center)[
  #block(
    fill: rgb("#6bb091"),
    inset: 8pt,
    radius: 5pt,
  )[
    $
    Y_(b,g,c) " " &~   "NegBinom"(mu_(b,g), phi.alt_g) \
    mu_(b,g) &= s_b dot q_g \
    arrow(1) dot arrow(s) &= 1
    $
  ]
]



We prove both statements below:

#pagebreak()
#block(
  fill: luma(235),
  inset: 8pt,
  radius: 5pt,
)[
  === The Constrained Model Is Identifiable
  Let $cal(M) := {pi_bold(theta) : bold(theta) in Theta}$ be a parameterized model where:
  - $bold(theta) := (arrow(s), arrow(q), arrow(phi.alt))$ with the constraint $arrow(1) dot arrow(s) = 1$.
  - $theta_(b,g) = (s_b, q_g, phi.alt_g)$.
  - $pi_theta_(b,g) := P_("NB")(n | s_b dot q_g, phi.alt_g)$.
  
  Let 
  $bold(theta) = (arrow(s), arrow(q), arrow(phi.alt))$ 
  and 
  $bold(theta') = (arrow(s'), arrow(q'), arrow(phi.alt'))$ with equality $pi_bold(theta) = pi_bold(theta')$ if and only if $pi_theta_(b,g) = pi_theta'_(b,g)$ for all $b,g$.
  
  Note that if two negative binomial distributions are equal, then both their realized location parameters $mu$ and their inverse overdispersion factors $phi.alt$ are equal. Since we are given $pi_bold(theta) = pi_bold(theta')$, it is immediate that $phi.alt_g = phi.alt'_g$ for all $g$. Thus, we must prove that equality of the location parameters $mu_(b,g) = mu'_(b,g) = s_b dot q_g =  s'_b dot q'_g$ implies $(s_b, q_g) = (s'_b, q'_g)$ for all $b,g$.
  
  The realization of the location parameter for all $b,g$ can be rewritten as the following outerproduct $bold(Mu) = arrow(s) dot arrow(q)^top$ with the $b$-th row and $g$-th column storing $mu_(b,g)$. Thus, equality of location parameters is equivalent to:
  $
  bold(Mu)
   &=
  bold(Mu') \
  arrow(s) dot arrow(q)^top
   &=
  arrow(s') dot arrow(q')^top \
  arrow(1) dot arrow(s) dot arrow(q)^top
   &=
  arrow(1) dot arrow(s') dot arrow(q')^top \
  arrow(q)^top
   &=
  arrow(q')^top \
  $
  Immediately implying $(arrow(s), arrow(q)) = (arrow(s'), arrow(q'))$.
]

#block(
  fill: luma(235),
  inset: 8pt,
  radius: 5pt,
)[
  === The Constrained Model Preserves Relative Differences
  Our assumed "true model" for the data generating process $cal(M) := {pi_bold(theta) : bold(theta) in Theta}$ and its constrained counterpart  $cal(M') := {pi_bold(theta') : bold(theta') in Theta'}$ are essentially equivalent. In order for relative differences to be preserved, then for $pi_bold(theta) in cal(M), pi_bold(theta') in cal(M')$ and $k in bb(R)$ :
  
  $
  pi_bold(theta) = pi_bold(theta') 
   ==> 
  eta_g = k q'_g quad quad forall g in {1, ..., n_g}
  $
  
  Let $pi_bold(theta) in cal(M), pi_bold(theta') in cal(M')$ be arbitrary with:
  - $bold(theta) := (arrow(zeta), arrow(eta), arrow(phi.alt))$ and $bold(theta') := (arrow(s'), arrow(q'), arrow(phi.alt'))$ with the constraint $arrow(1) dot arrow(s) = 1$. 
  - $pi_theta_(b,g) := P_("NB")(n | zeta_b dot eta_g, phi.alt_g)$ and $pi_theta'_(b,g) := P_("NB")(n | s'_b dot q'_g, phi.alt'_g)$.
  - $b in {1, dots, n_b}$ and $g in {1, dots, n_g}$.
  Then $pi_bold(theta) = pi_bold(theta')$ implies:
  
  $
  pi_bold(theta) = pi_bold(theta')
   &==>
  P_("NB")(n | zeta_b dot eta_g, phi.alt_g) = P_("NB")(n | s'_b dot q'_g, phi.alt'_g) quad forall b,g \
   &==>
  zeta_b dot eta_g = s'_b dot q'_g quad forall b,g \
   &==>
  arrow(zeta) dot arrow(eta)^top = arrow(s) dot arrow(q)^top \
   &==>
  arrow(1) dot arrow(zeta) dot arrow(eta)^top = arrow(1) dot arrow(s) dot arrow(q)^top \
   &==>
  arrow(1) dot arrow(zeta) dot arrow(eta)^top = arrow(q)^top \
   &==>
  arrow(eta)^top = frac(1, arrow(1) dot arrow(zeta)) arrow(q)^top 
  $
  
  Therefore, there exists some constant $k = frac(1, arrow(1) dot arrow(zeta))$ such that $pi_bold(theta) = pi_bold(theta') ==> eta_g = k q'_g quad forall g in {1, ..., n_g}$.
]

== Modelling Multiple Cell Populations From a Single Tissue

The natural extension to multiple cell populations is introduced by adding another index variable $p$ to group observations by their annotated cell-type:

#align(center)[
  #block(
    fill: rgb("#6bb091"),
    inset: 8pt,
    radius: 5pt,
  )[
    $
    Y_(b,p,g,c) " " &~  "NegBinom"(mu_(b,p,g), phi.alt_g) \
    mu_(b,p,g) &= s_b dot q_(p,g) \
    arrow(1) dot arrow(s) &= 1
    $
  ]
]

Note however that there is potential for this model to be unidentifiable if there exists a batch that by chance contains no cell populations shared with other batches. E.g, if you seek to compare expression for two populations that were processed separately, you cannot infer what differences are due to batch-effects and what differences are due to biology. Fortunately, this scenario is unlikely all samples are derived from the same tissue--- it is expected that the same cell populations will repeatedly show up.

== Modelling Multiple Cell Populations of a Single Tissue Across Donors

The extension to multiple donors is also introduced by adding another index variable $d$ to group observations by the donor they were sampled from:

#align(center)[
  #block(
    fill: rgb("#C17985"),
    inset: 8pt,
    radius: 5pt,
  )[
    $
    Y_(b,d,p,g,c) " " &~  "NegBinom"(mu_(b,d,p,g), phi.alt_g) \
    mu_(b,d,p,g) &= s_b dot q_(d,p,g) \
    arrow(1) dot arrow(s) &= 1
    $
  ]
]

But an issue stems from the fact that individuals often differ in their biology in some way. Differences in genetics and environmental exposures inevitably lead to changes in expression even when looking at the same cell population. I.e, for two donors $d_1, d_2$ and a cell population $p$: $eta_(d_1, p, g)$ is not guarenteed to be equal to $eta_(d_2, p, g)$ for all genes. This complicates inference as samples from different donors are also often processsed in different batches, so we run into the same issue described above when two populations are processed separately--- we cannot distinguish between biology and batch-effect.

We need insight on either the size factors or the gene expression to make our model identifiable again. The ideal method would be to add "molecular spikes"(artificial RNAs with known sequence and fixed concentration) into each batch prior to sequencing(anecdotally the major cause behind batch-effects), and estimate a common intercept $alpha$ that is fixed across batches for the number of reads $S_b$ mapping to the spike-in sequence:

#align(center)[
  #block(
    fill: rgb("#6bb091"),
    inset: 8pt,
    radius: 5pt,
  )[
    $
    Y_(b,d,p,g,c) " " &~  "NegBinom"(mu_(b,d,p,g), phi.alt_g) \
    mu_(b,d,p,g) &= s_b dot q_(d,p,g) \
    arrow(1) dot arrow(s) &= 1
    \
    \
    S_b " " &~  "Poisson"(mu_b) \
    mu_b &= s_b dot alpha
    $
  ]
]

The authors behind Tabula Sapiens however did not use molecular spikes. So we must either give up or make assumptions about gene expression. Fortunately, we have good reason to make assumptions about gene expression.

#pagebreak()

Despite differences in biology, it is plausible that the expression of a gene $g$ is more similar among donors _within_ a particular cell population $p$ compared to _between_ populations. 

We adopt the assumption that for a fixed cell population $p$ and gene $g$, expression quantities are log-normally distributed across donors with parameters $(alpha_(p,g), lambda_(p,g))$. It should be noted that the specific choice of distribution here is somewhat arbitrary; any distribution that provides control over the mean and variance of expression quantities could be used.

$
q_(d,p,g) " " &~ "Log-Normal"(alpha_(p,g), lambda_(p,g) tau_p) \
lambda_(p,g) " " &~ " Half-Cauchy"^+(0, 1)
$

The similarity of a cell population's gene expression across donors is represented by $lambda_(p,g)$--- approximately zero for genes whose variation is primarily due to batch-effects. Following @carvalho_handling_2009, we model the shrinkage of this variance term by including a half-Cauchy prior:

$

$

Where $tau_p$ represents the global shrinkage parameter.

== Bells and Whistles

The bulk of the model has been designed above step-by-step, gradually increasing the complexity and working around issues that pop up. The two final additions are 1) a log-normal prior on the inverse overdispersion factors $phi.alt_g$ to partially pool estimates and 2) an exponential prior on the per-population global shrinkage parameter $tau_p$ with $c^* = 1$ by default:

#align(center)[
  #block(
    fill: rgb("#6bb091"),
    inset: 8pt,
    radius: 5pt,
  )[
    $
    Y_(b,d,p,g,c) " " &~  "NegBinom"(mu_(b,d,p,g), phi.alt_g) \
    mu_(b,d,p,g) &= s_b dot q_(d,p,g) \
    arrow(1) dot arrow(s) &= 1 \
    \
    q_(d,p,g) " " &~  "Log-Normal"(alpha_(p,g), lambda_(p,g)) \
    lambda_(p,g) " " &~  " Half-Cauchy" ^+ (0, tau_p) \
    tau_p " " &~  "Exponential"(c^*) \
    \
    phi.alt_g " " &~  "Log-Normal"(psi, sigma)
    $
  ]
]

#pagebreak()
= Validation



== Modelling a Single Population Across 2 Donors

Simulation parameters:
- Number of genes : 1000
- Number of genes with nonzero differences : 10
- Number of cells per-donor : 200

#figure(
  image(
    "../src/modelling/simulations/single-population/plot_contrasts.png",
    width: 95%),
  caption: [95% credibility intervals for the ratio of relative expression quantities between donors `1` and `0`. Coloured dots represent the ground truth ratios(black) and size factor-adjusted sample mean ratios(grey). Only the first 100 transcripts are shown here.]
)

#pagebreak()
== Modelling Multiple Populations Across 2 Donors

Simulation parameters:
- Number of populations: 3
- Number of genes : 500
- Number of genes with nonzero differences : 10
- Number of cells per-population : 200

#figure(
  image(
    "../src/modelling/simulations/multiple-populations/plot_contrasts.png",
    width: 95%),
  caption: [95% credibility intervals for the ratio of relative expression quantities between populations `2` and `1`. Coloured dots represent the ground truth ratios(black) and size factor-adjusted sample mean ratios(grey). Only the first 100 transcripts are shown here.]
)

= An Outlook on Previous Methods

There is no common method for the joint modelling of scRNA-seq data for differential expression analysis as I've described above. The closest that I know of is ZINBMM @li_zinbmm_2023, which is essentially a finite mixture model with a zero-inflated negative binomial distributional assumption. Disregarding the zero-inflation assumption(see @svensson_droplet_2020 for some commentary if you're interested), I do like the idea of fitting a finite mixture model on scRNA-seq data. Unfortunately, 1) there are often continuous trajectories that can mess up the identification of discrete cell populations, and 2) most people are already used to clustering being separate from differential expression analysis.

After clustering, the most common procedure for differential expression analysis is to pseudobulk(sum counts per-gene, per-batch, per-condition, etc) cell types of interest and compute however many pairwise comparisons with existing bulk-RNAseq methods(edgeR, DESeq2, etc). These methods use heuristics to calculate size factors, so I would be interested in seeing how well they match up to my proposed method.

I also realized a few things while writing this, it would be relatively easy to extend the above model and handle arbitrary covariates. The only reason donor($d$) and population($p$) are the main focus here is because the original goal was to infer cell type-specific expression quantities without molecular spikes. I think it would be interesting to try and write a general framework for jointly modelling RNAseq data, but the major appeal I see is that I'm not sure if pseudo-bulking is the right choice particularly for scRNA-seq(it's better to just treat them as repeated measurements/i.i.d observations rather than summing because I think the number of cells collected will bias you in some way, last I checked people are still summing) so that might be a knowledge gap to build a tool for or at least correct people.

On a final note, I find the idea of placing a prior on a subset of genes that are more likely to be conserved interesting. For example, if we know _a priori_ that a certain set of genes are known to not differ between two or more groups(whether the groups be donors, populations, conditions, or otherwise), then this is useful information for estimating relative size factors.


#bibliography("refs.bib", full : true)