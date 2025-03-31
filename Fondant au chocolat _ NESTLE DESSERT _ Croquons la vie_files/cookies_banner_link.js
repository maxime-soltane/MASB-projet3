/**
 * @file
 * Javascript file for brand colors.
 */

(function ($, Drupal, drupalSettings) {

  Drupal.behaviors.clvCookieBannerLink = {
    attach: function (context, settings) {
      // Open cookie banner.
      $('.cookie-banner-link').on('click', function () {
        window.evidon.notice.showOptions(null);
        return false;
      });
    }
  };

})(jQuery, Drupal, drupalSettings);
