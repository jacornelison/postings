/**
 * Toggles the showhide_hidden style class. Utilizing styles in showhide.css,
 * the contents can be gracefully shown and hidden using CSS3 transformations.
 * 
 * showhide(selector)
 * 
 *     selector    A CSS-style selector which specifies the element to operate
 *                 upon.
 */
function showhide(selector)
{
    var el = document.querySelector(selector);
    el.classList.toggle("showhide_hidden");
}
