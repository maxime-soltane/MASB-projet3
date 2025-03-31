(function ($) {
  "use strict";

  Drupal.behaviors.bazaarvoiceInlineProductRatings = {
    attach: function (context, settings) {
      if(drupalSettings.bazaarvoiceProducts.apiVersion === 'bv') {
        return;
      }
      var product_ids = drupalSettings.bazaarvoiceProducts.productIds;
      $BV.ui('rr', 'inline_ratings', {
          productIds: product_ids,
          containerPrefix: 'BVRRInlineRating'
        }
      );
    }
  }
})(jQuery);