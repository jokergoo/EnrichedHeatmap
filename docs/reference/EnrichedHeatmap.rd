<!DOCTYPE html>
<!-- Generated by pkgdown: do not edit by hand --><html lang="en"><head><meta http-equiv="Content-Type" content="text/html; charset=UTF-8"><meta charset="utf-8"><meta http-equiv="X-UA-Compatible" content="IE=edge"><meta name="viewport" content="width=device-width, initial-scale=1.0"><title>Constructor Method for the Enriched Heatmap — EnrichedHeatmap • EnrichedHeatmap</title><!-- jquery --><script src="https://cdnjs.cloudflare.com/ajax/libs/jquery/3.4.1/jquery.min.js" integrity="sha256-CSXorXvZcTkaix6Yvo6HppcZGetbYMGWSFlBw8HfCJo=" crossorigin="anonymous"></script><!-- Bootstrap --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/css/bootstrap.min.css" integrity="sha256-bZLfwXAP04zRMK2BjiO8iu9pf4FbLqX6zitd+tIvLhE=" crossorigin="anonymous"><script src="https://cdnjs.cloudflare.com/ajax/libs/twitter-bootstrap/3.4.1/js/bootstrap.min.js" integrity="sha256-nuL8/2cJ5NDSSwnKD8VqreErSWHtnEP9E7AySL+1ev4=" crossorigin="anonymous"></script><!-- bootstrap-toc --><link rel="stylesheet" href="../bootstrap-toc.css"><script src="../bootstrap-toc.js"></script><!-- Font Awesome icons --><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/all.min.css" integrity="sha256-mmgLkCYLUQbXn0B1SRqzHar6dCnv9oZFPEC1g1cwlkk=" crossorigin="anonymous"><link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/5.12.1/css/v4-shims.min.css" integrity="sha256-wZjR52fzng1pJHwx4aV2AO3yyTOXrcDW7jBpJtTwVxw=" crossorigin="anonymous"><!-- clipboard.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/clipboard.js/2.0.6/clipboard.min.js" integrity="sha256-inc5kl9MA1hkeYUt+EC3BhlIgyp/2jDIyBLS6k3UxPI=" crossorigin="anonymous"></script><!-- headroom.js --><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/headroom.min.js" integrity="sha256-AsUX4SJE1+yuDu5+mAVzJbuYNPHj/WroHuZ8Ir/CkE0=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/headroom/0.11.0/jQuery.headroom.min.js" integrity="sha256-ZX/yNShbjqsohH1k95liqY9Gd8uOiE1S4vZc+9KQ1K4=" crossorigin="anonymous"></script><!-- pkgdown --><link href="../pkgdown.css" rel="stylesheet"><script src="../pkgdown.js"></script><meta property="og:title" content="Constructor Method for the Enriched Heatmap — EnrichedHeatmap"><meta property="og:description" content="Constructor Method for the Enriched Heatmap"><!-- mathjax --><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js" integrity="sha256-nvJJv9wWKEm88qvoQl9ekL2J+k/RWIsaSScxxlsrv8k=" crossorigin="anonymous"></script><script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/config/TeX-AMS-MML_HTMLorMML.js" integrity="sha256-84DKXVJXs0/F8OTMzX4UR909+jtl4G7SPypPavF+GfA=" crossorigin="anonymous"></script><!--[if lt IE 9]>
<script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
<script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
<![endif]--></head><body data-spy="scroll" data-target="#toc">
    

    <div class="container template-reference-topic">
      <header><div class="navbar navbar-default navbar-fixed-top" role="navigation">
  <div class="container">
    <div class="navbar-header">
      <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false">
        <span class="sr-only">Toggle navigation</span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
        <span class="icon-bar"></span>
      </button>
      <span class="navbar-brand">
        <a class="navbar-link" href="../index.html">EnrichedHeatmap</a>
        <span class="version label label-default" data-toggle="tooltip" data-placement="bottom" title="">1.33.1</span>
      </span>
    </div>

    <div id="navbar" class="navbar-collapse collapse">
      <ul class="nav navbar-nav"><li class="dropdown">
  <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" data-bs-toggle="dropdown" aria-expanded="false">
    Articles
     
    <span class="caret"></span>
  </a>
  <ul class="dropdown-menu" role="menu"><li class="dropdown-header">Vignettes</li>
    <li>
      <a href="../articles/EnrichedHeatmap_intro.html">Make Enriched Heatmaps</a>
    </li>
    <li>
      <a href="../articles/visualize_categorical_signals.html">Visualize Categorical Signals</a>
    </li>
    <li>
      <a href="../articles/row_ordering.html">Compare row ordering methods</a>
    </li>
    <li>
      <a href="../articles/roadmap.html">Visualize Comprehensive Associations in Roadmap dataset</a>
    </li>
    <li class="divider">
    <li class="dropdown-header">More applications</li>
  </ul></li>
<li>
  <a href="../reference/index.html">Reference</a>
</li>
      </ul><ul class="nav navbar-nav navbar-right"><li>
  <a href="https://github.com/jokergoo/EnrichedHeatmap/" class="external-link">
    <span class="fab fa-github fa-lg"></span>
     
  </a>
</li>
      </ul></div><!--/.nav-collapse -->
  </div><!--/.container -->
</div><!--/.navbar -->

      

      </header><div class="row">
  <div class="col-md-9 contents">
    <div class="page-header">
    <h1>Constructor Method for the Enriched Heatmap</h1>
    
    <div class="hidden name"><code>EnrichedHeatmap.rd</code></div>
    </div>

    <div class="ref-description">
    <p>Constructor Method for the Enriched Heatmap</p>
    </div>

    <div id="ref-usage">
    <div class="sourceCode"><pre class="sourceCode r"><code><span><span class="fu">EnrichedHeatmap</span><span class="op">(</span><span class="va">mat</span>,</span>
<span>    <span class="va">col</span>,</span>
<span>    top_annotation <span class="op">=</span> <span class="fu">HeatmapAnnotation</span><span class="op">(</span>enriched <span class="op">=</span> <span class="fu"><a href="anno_enriched.html">anno_enriched</a></span><span class="op">(</span><span class="op">)</span><span class="op">)</span>,</span>
<span>    row_order <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/order.html" class="external-link">order</a></span><span class="op">(</span><span class="fu"><a href="enriched_score.html">enriched_score</a></span><span class="op">(</span><span class="va">mat</span><span class="op">)</span>, decreasing <span class="op">=</span> <span class="cn">TRUE</span><span class="op">)</span>,</span>
<span>    pos_line <span class="op">=</span> <span class="cn">TRUE</span>,</span>
<span>    pos_line_gp <span class="op">=</span> <span class="fu">gpar</span><span class="op">(</span>lty <span class="op">=</span> <span class="fl">2</span><span class="op">)</span>,</span>
<span>    axis_name <span class="op">=</span> <span class="cn">NULL</span>,</span>
<span>    axis_name_rot <span class="op">=</span> <span class="fl">0</span>,</span>
<span>    axis_name_gp <span class="op">=</span> <span class="fu">gpar</span><span class="op">(</span>fontsize <span class="op">=</span> <span class="fl">10</span><span class="op">)</span>,</span>
<span>    border <span class="op">=</span> <span class="cn">TRUE</span>,</span>
<span>    cluster_rows <span class="op">=</span> <span class="cn">FALSE</span>,</span>
<span>    row_dend_reorder <span class="op">=</span> <span class="op">-</span><span class="fu"><a href="enriched_score.html">enriched_score</a></span><span class="op">(</span><span class="va">mat</span><span class="op">)</span>,</span>
<span>    show_row_dend <span class="op">=</span> <span class="cn">FALSE</span>,</span>
<span>    show_row_names <span class="op">=</span> <span class="cn">FALSE</span>,</span>
<span>    heatmap_legend_param <span class="op">=</span> <span class="fu"><a href="https://rdrr.io/r/base/list.html" class="external-link">list</a></span><span class="op">(</span><span class="op">)</span>,</span>
<span>    <span class="va">...</span><span class="op">)</span></span></code></pre></div>
    </div>

    <div id="arguments">
    <h2>Arguments</h2>
    <dl><dt>mat</dt>
<dd><p>A matrix which is returned by <code><a href="normalizeToMatrix.html">normalizeToMatrix</a></code>.</p></dd>

  <dt>col</dt>
<dd><p>Color settings. If the signals are categorical, color should be a vector with category levels as names.</p></dd>

  <dt>top_annotation</dt>
<dd><p>A special annotation which is always put on top of the enriched heatmap and is constructed by <code><a href="anno_enriched.html">anno_enriched</a></code>.</p></dd>

  <dt>row_order</dt>
<dd><p>Row order. Default rows are ordered by enriched scores calculated from <code><a href="enriched_score.html">enriched_score</a></code>.</p></dd>

  <dt>pos_line</dt>
<dd><p>Whether draw vertical lines which represent the positions of <code>target</code>?</p></dd>

  <dt>pos_line_gp</dt>
<dd><p>Graphic parameters for the position lines.</p></dd>

  <dt>axis_name</dt>
<dd><p>Names for axis which is below the heatmap. If the targets are single points, <code>axis_name</code> is a vector of length three which corresponds to upstream, target itself and downstream. If the targets are regions with width larger than 1, <code>axis_name</code> should be a vector of length four which  corresponds to upstream, start of targets, end of targets and downstream.</p></dd>

  <dt>axis_name_rot</dt>
<dd><p>Rotation for axis names.</p></dd>

  <dt>axis_name_gp</dt>
<dd><p>Graphic parameters for axis names.</p></dd>

  <dt>border</dt>
<dd><p>Whether show the border of the heatmap?</p></dd>

  <dt>cluster_rows</dt>
<dd><p>Clustering on rows are turned off by default.</p></dd>

  <dt>show_row_dend</dt>
<dd><p>Whether show dendrograms on rows if hierarchical clustering is applied on rows?</p></dd>

  <dt>row_dend_reorder</dt>
<dd><p>Weight for reordering the row dendrogram. It is reordered by enriched scores by default.</p></dd>

  <dt>show_row_names</dt>
<dd><p>Whether show row names?</p></dd>

  <dt>heatmap_legend_param</dt>
<dd><p>A list of settings for heatmap legends. <code>at</code> and <code>labels</code> can not be set here.</p></dd>

  <dt>...</dt>
<dd><p>Other arguments passed to <code><a href="https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html" class="external-link">Heatmap</a></code>.</p></dd>


</dl></div>
    <div id="details">
    <h2>Details</h2>
    <p>The enriched heatmap is essentially a normal heatmap but with several special settings. Following parameters are 
set with pre-defined values:</p>
<dl><dt><code>cluster_columns</code></dt>
<dd><p>enforced to be <code>FALSE</code></p></dd>

  <dt><code>show_column_names</code></dt>
<dd><p>enforced to be <code>FALSE</code></p></dd>

  <dt><code>bottom_annotation</code></dt>
<dd><p>enforced to be <code>NULL</code></p></dd>


</dl><p><code>EnrichedHeatmap</code> calls <code><a href="https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html" class="external-link">Heatmap</a></code>, thus, most of the 
arguments in <code><a href="https://rdrr.io/pkg/ComplexHeatmap/man/Heatmap.html" class="external-link">Heatmap</a></code> are usable in <code>EnrichedHeatmap</code> such as
to apply clustering on rows, or to split rows by a data frame or k-means clustering. Users can also 
add more than one heatmaps by <code>+</code> operator. Enriched heatmaps and normal heatmaps can be
concatenated mixed.</p>
<p>For detailed demonstration, please go to the vignette.</p>
    </div>
    <div id="value">
    <h2>Value</h2>
    

<p>A <code>Heatmap-class</code> object.</p>
    </div>
    <div id="author">
    <h2>Author</h2>
    <p>Zuguang Gu &lt;z.gu@dkfz.de&gt;</p>
    </div>

    <div id="ref-examples">
    <h2>Examples</h2>
    <div class="sourceCode"><pre class="sourceCode r"><code><span class="r-in"><span><span class="fu"><a href="https://rdrr.io/r/base/load.html" class="external-link">load</a></span><span class="op">(</span><span class="fu"><a href="https://rdrr.io/r/base/system.file.html" class="external-link">system.file</a></span><span class="op">(</span><span class="st">"extdata"</span>, <span class="st">"chr21_test_data.RData"</span>, package <span class="op">=</span> <span class="st">"EnrichedHeatmap"</span><span class="op">)</span><span class="op">)</span></span></span>
<span class="r-in"><span><span class="va">mat3</span> <span class="op">=</span> <span class="fu"><a href="normalizeToMatrix.html">normalizeToMatrix</a></span><span class="op">(</span><span class="va">meth</span>, <span class="va">cgi</span>, value_column <span class="op">=</span> <span class="st">"meth"</span>, mean_mode <span class="op">=</span> <span class="st">"absolute"</span>,</span></span>
<span class="r-in"><span>    extend <span class="op">=</span> <span class="fl">5000</span>, w <span class="op">=</span> <span class="fl">50</span>, smooth <span class="op">=</span> <span class="cn">TRUE</span><span class="op">)</span></span></span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> All signal values are within [0, 1], so we assume it is methylation</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> signal. Automatically set limit [0, 1] to the smoothed values. If this</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> is not the case, set argument `limit = NA` in the function to remove</span>
<span class="r-msg co"><span class="r-pr">#&gt;</span> the limits. Set `verbose = FALSE` to turn off this message.</span>
<span class="r-in"><span><span class="fu">EnrichedHeatmap</span><span class="op">(</span><span class="va">mat3</span>, name <span class="op">=</span> <span class="st">"methylation"</span>, column_title <span class="op">=</span> <span class="st">"methylation near CGI"</span><span class="op">)</span></span></span>
<span class="r-plt img"><img src="EnrichedHeatmap-1.png" alt="" width="700" height="433"></span>
<span class="r-in"><span><span class="fu">EnrichedHeatmap</span><span class="op">(</span><span class="va">mat3</span>, name <span class="op">=</span> <span class="st">"meth1"</span><span class="op">)</span> <span class="op">+</span> <span class="fu">EnrichedHeatmap</span><span class="op">(</span><span class="va">mat3</span>, name <span class="op">=</span> <span class="st">"meth2"</span><span class="op">)</span></span></span>
<span class="r-plt img"><img src="EnrichedHeatmap-2.png" alt="" width="700" height="433"></span>
<span class="r-in"><span><span class="co"># for more examples, please go to the vignette</span></span></span>
</code></pre></div>
    </div>
  </div>
  <div class="col-md-3 hidden-xs hidden-sm" id="pkgdown-sidebar">
    <nav id="toc" data-toggle="toc" class="sticky-top"><h2 data-toc-skip>Contents</h2>
    </nav></div>
</div>


      <footer><div class="copyright">
  <p></p><p>Developed by Zuguang Gu.</p>
</div>

<div class="pkgdown">
  <p></p><p>Site built with <a href="https://pkgdown.r-lib.org/" class="external-link">pkgdown</a> 2.0.7.</p>
</div>

      </footer></div>

  


  <style>nav[data-toggle='toc'] .nav .nav {display: block;}</style></body></html>

