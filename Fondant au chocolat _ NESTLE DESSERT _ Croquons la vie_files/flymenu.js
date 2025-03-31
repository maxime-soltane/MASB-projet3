/**
 * @file
 * Handles Flymenu integration.
 */

(function ($, Drupal, drupalSettings, once) {

  /**
   * @type {{attach: Drupal.behaviors.showFlyMenuCart.attach}}
   */
  Drupal.behaviors.showFlyMenuCart = {
    attach: function (context, settings) {
      // Initialisation FlyMenu.
      if (typeof flyMenuInit !== 'undefined') {
        flyMenuInit({
          site_token: drupalSettings.FlyMenuGlobalInfo.site_token,
          flymenuWidgetPosition: drupalSettings.FlyMenuGlobalInfo.widget_position,
          session_lifetime: drupalSettings.FlyMenuGlobalInfo.session_lifetime,
          show_recipe_page: drupalSettings.FlyMenuGlobalInfo.show_recipe_page,
          showOrderButtonOnRecipeList: drupalSettings.FlyMenuGlobalInfo.show_order_button
        });

        $(once('processed', ".flymenubutton")).click(function () {
          Drupal.behaviors.showFlyMenuCart.addRecipe();
        });
      }
    },
    addRecipe: function () {
      if (typeof addRecipes !== 'undefined') {
        addRecipes({
          recipes: drupalSettings.FlyMenuRecipeInfo
        });
      }
    },
  };
})(jQuery, Drupal, drupalSettings, once);
