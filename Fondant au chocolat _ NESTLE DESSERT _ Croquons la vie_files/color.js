/**
 * @file
 * Javascript file for brand colors.
 */

(function ($, Drupal, drupalSettings, once) {

  Drupal.behaviors.clvBrandColor = {
    attach: function (context, settings) {
      $(once('processed', 'head')).each(function () {
        var css = '<style>:root {';
        if (drupalSettings.colorBackgroundLogo !== null) {
          css += '--color-background-logo: ' + drupalSettings.colorBackgroundLogo + ';';
        }
        if (drupalSettings.colorBrand !== null) {
          css += '--color-brand: ' + drupalSettings.colorBrand + ';';
        }
        if (drupalSettings.colorSecondary !== null) {
          css += '--color-secondary: ' + drupalSettings.colorSecondary + ';';
        }
        if (drupalSettings.colorBrandText !== null) {
          css += '--color-brand-text: ' + drupalSettings.colorBrandText + ';';
        }
        if (drupalSettings.colorSecondaryText !== null) {
          css += '--color-secondary-text: ' + drupalSettings.colorSecondaryText + ';';
        }
        css += '}</style>';
        $(this).append(css);
      });
    }
  };

})(jQuery, Drupal, drupalSettings, once);
