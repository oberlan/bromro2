# Instances

Instances used in the papers:

* `bromro`: [The multiperiod two-dimensional non-guillotine cutting stock problem with usable leftovers](https://onlinelibrary.wiley.com/doi/abs/10.1111/itor.12648)
* `bromro2`: [A forward-looking matheuristic approach for the multi-period two-dimensional non-guillotine cutting stock problem with usable leftovers](https://www.ime.usp.br/~egbirgin/)


## Exemple of instance:

```
4

3	8	10	9	
11	10	19	
8	18	10	
19	19	19	18	7	7	9	18	
13	13	13	7	6	6	6	18	
17	17	11	11	10	10	10	18	18	18	
10	10	5	5	16	16	16	7	7	7	
7	7	7	20	20	20	14	14	14	
19	19	19	19	19	19	7	7	7	


2
7	11	
6	5	

3	3	3	3	
59	79	33	
86	71	40	
1	1	1	
69	85	62	
54	43	69	
1	1	1	
59	59	96	
32	37	43	
1	1	1	
46	39	39	
54	64	64	
1	1	1	
```

The first value represent the number of periods ($P$). The next line contains $P$ values indicating the number of items ($n_s$, $s = 1, \dots, P$). The next $2 \times P$ lines contains the dimensions of each item for each period, following the pattern:

$w_{11}~~w_{12}~~\dots~~w_{1n_1}$
$h_{11}~~h_{12}~~\dots~~h_{1n_1}$
$w_{21}~~w_{22}~~\dots~~w_{2n_2}$
$h_{21}~~h_{22}~~\dots~~h_{2n_2}$
$\vdots$
$w_{P1}~~w_{P2}~~\dots~~w_{Pn_P}$
$h_{P1}~~h_{P2}~~\dots~~h_{Pn_P}$

The next line contains the amount of catalog itens ($d$). The next two lines contains the dimensions of each item in the catalog, following the pattern:
$\bar{w}_1~~\bar{w}_1~~\dots~~\bar{w}_d$
$\bar{h}_1~~\bar{h}_1~~\dots~~\bar{h}_d$

Next, the amount of available objects ($m_s$, $s = 1, \dots, P$) in each period is shown, following the pattern:

$W_{11}~~W_{12}~~\dots~~W_{1m_1}$
$H_{11}~~H_{12}~~\dots~~H_{1m_1}$
$c_{11}~~c_{12}~~\dots~~c_{1m_1}$
$W_{21}~~W_{22}~~\dots~~W_{2m_2}$
$H_{21}~~H_{22}~~\dots~~H_{2m_2}$
$c_{21}~~c_{22}~~\dots~~c_{2m_2}$
$\vdots$
$W_{P1}~~W_{P2}~~\dots~~W_{Pm_P}$
$H_{P1}~~H_{P2}~~\dots~~H_{Pm_P}$
$c_{P1}~~c_{P2}~~\dots~~c_{Pm_P}$

In summary, the instance follows the pattern:

$P$

$n_1~~n_2~~\dots~~n_P$
$w_{11}~~w_{12}~~\dots~~w_{1n_1}$
$h_{11}~~h_{12}~~\dots~~h_{1n_1}$
$w_{21}~~w_{22}~~\dots~~w_{2n_2}$
$h_{21}~~h_{22}~~\dots~~h_{2n_2}$
$\vdots$
$w_{P1}~~w_{P2}~~\dots~~w_{Pn_P}$
$h_{P1}~~h_{P2}~~\dots~~h_{Pn_P}$

$d$
$\bar{w}_1~~\bar{w}_1~~\dots~~\bar{w}_d$
$\bar{h}_1~~\bar{h}_1~~\dots~~\bar{h}_d$

$m_1~~m_2~~\dots~~m_P$
$W_{11}~~W_{12}~~\dots~~W_{1m_1}$
$H_{11}~~H_{12}~~\dots~~H_{1m_1}$
$c_{11}~~c_{12}~~\dots~~c_{1m_1}$
$W_{21}~~W_{22}~~\dots~~W_{2m_2}$
$H_{21}~~H_{22}~~\dots~~H_{2m_2}$
$c_{21}~~c_{22}~~\dots~~c_{2m_2}$
$\vdots$
$W_{P1}~~W_{P2}~~\dots~~W_{Pm_P}$
$H_{P1}~~H_{P2}~~\dots~~H_{Pm_P}$
$c_{P1}~~c_{P2}~~\dots~~c_{Pm_P}$