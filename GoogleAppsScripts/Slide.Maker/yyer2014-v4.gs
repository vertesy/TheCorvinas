/*  Made by: yyer2014@gmail.com (Alexander Y.)

  Main conversion function.
  Top level list items work as new slide titles. Slide layout is TITLE_AND_BODY
  Top level headings make slides with TITLE layout.
  Subheadings (levels 1 and 2) produce SECTION_HEADER slide layouts.
  Inline images are placed on the corresponding slide and centered.

  The source document has SOURCE_ID and the target presentation has TARGET_ID Google identifiers.
*/
function docs2slides() {
  const SOURCE_ID = '  DOCUMENT ID  ',
      TARGET_ID = '  DOCUMENT ID  '

  // Open target presentation and remember old slides count (they should be removed)
  var slides = SlidesApp.openById(TARGET_ID),
      slidesCount = slides.getSlides().length;

  var body = DocumentApp.openById(SOURCE_ID).getBody(),
      data = {
        slides: slides,
        body: [],
        images: []
      };
  var headings = {};
  headings[DocumentApp.ParagraphHeading.HEADING1] = SlidesApp.PredefinedLayout.TITLE;
  headings[DocumentApp.ParagraphHeading.HEADING2] = SlidesApp.PredefinedLayout.SECTION_HEADER;
  headings[DocumentApp.ParagraphHeading.HEADING3] = SlidesApp.PredefinedLayout.SECTION_HEADER;

  // Loop through all source body child elements to find text paragraphs or list items
  // As a result data object will include h1 member - for a slide title, and body member - for a slide rich text content.
  // Ready data object is the only argument to call addSlide function for a new slide creation
  for (var ci = 0; ci < body.getNumChildren(); ci++) {
    var child = body.getChild(ci);
    switch(child.getType()) {
      case DocumentApp.ElementType.PARAGRAPH:
        var p = child.asParagraph();
        // Paragraph may have a heading type
        if (headings[p.getHeading()] != undefined) {
          addSlide(data);
          data.layout = headings[p.getHeading()];
          data.h1 = p.getText();  // This is string type
          data.body = [];
          data.images = [];
        } else {  // or normal type (no heading)
          data.body.push(p.editAsText());  // This is rich Text type (NOT string!)
        }
        appendImagesFrom(p, data);
        break;
      case DocumentApp.ElementType.LIST_ITEM:
        var li = child.asListItem();
        // List item of the top level is considered as a slide title
        if (li.getNestingLevel() == 0) {
          addSlide(data);
          data.layout = SlidesApp.PredefinedLayout.TITLE_AND_BODY;
          data.h1 = li.getText();  // string type
          data.body = [];
          data.images = [];
        } else {
          data.body.push(li.editAsText());
        }
        appendImagesFrom(li, data);
    }
  }
  addSlide(data);

  // Remove old slides (if new slides exist)
  while (slides.getSlides().length > slidesCount && slidesCount > 0) {
    slides.getSlides()[0].remove();
    slidesCount--;
  }

}


function appendImagesFrom(parentElement, toData) {
  for (var i = 0; i < parentElement.getNumChildren(); i++) {
    var c = parentElement.getChild(i);
    if (c.getType() == DocumentApp.ElementType.INLINE_IMAGE) {
      toData.images.push(c);
    }
  }
}


/* This function takes data object which has members:
h1 - slide title,
body - slide rich text content
layout - desired slide layout
images - inline images found for the slide
slides - presentation itself

It appends a slide with the predefined layout, fills title and content, applies styles to content.
The styles are discovered from data.body array elements, because they have rich Text type
*/
function addSlide(data) {
  if (data.h1 == undefined) return;
  var slide = data.slides.appendSlide(data.layout);
  slide.getShapes().forEach(function(shape) {
    switch(shape.getPlaceholderType()) {
      case SlidesApp.PlaceholderType.TITLE:
      case SlidesApp.PlaceholderType.CENTERED_TITLE:
        shape.getText().setText(data.h1);  // fill in a slide title
        break;
      case SlidesApp.PlaceholderType.BODY:
      case SlidesApp.PlaceholderType.SUBTITLE:
        var textRng = shape.getText();
        data.body.forEach(
          function(text) {
            var t = text.getText();  // string type
            var p = textRng.appendParagraph(t);
            // This is important call to discover styles (rich text formatting) inside text elements
            var styles = discoverStyles(text);
            // Once discovered the styles should be applied to text fragments accordingly
            styles[0].forEach(function(s) {
              if (s.end == undefined) s.end = t.length;
              p.getRange().getRange(s.start, s.end).getTextStyle().setBold(true);
            });
            styles[1].forEach(function(s) {
              if (s.end == undefined) s.end = t.length;
              p.getRange().getRange(s.start, s.end).getTextStyle().setItalic(true);
            });
            styles[2].forEach(function(s) {
              if (s.end == undefined) s.end = t.length;
              p.getRange().getRange(s.start, s.end).getTextStyle().setUnderline(true);
            });
            styles[3].forEach(function(s) {
              if (s.end == undefined) s.end = t.length;
              p.getRange().getRange(s.start, s.end).getTextStyle().setForegroundColor(s.value);
            });
            styles[4].forEach(function(s) {
              if (s.end == undefined) s.end = t.length;
              p.getRange().getRange(s.start, s.end).getTextStyle().setBackgroundColor(s.value);
            });
            styles[5].forEach(function(s) {
              if (s.end == undefined) s.end = t.length;
              p.getRange().getRange(s.start, s.end).getTextStyle().setLinkUrl(s.value);
            });
          }
        );
        break;
      default:
        Logger.log(shape.getPlaceholderType());
    }
  });
  data.images.forEach(function(image) {
    var blob = image.getBlob();
    slide.insertImage(blob).alignOnPage(SlidesApp.AlignmentPosition.CENTER);
  });
  data.h1 = undefined;
}


/* We are interested in 6 style features: bold, italic, underline, color, background color, hyperlink
All of them are considered as an ordered array with 6 elements. Such an array appears in different places:
  styles - final resulting structure to be returned
  v - the style vector for a particular position in text
  currentStyle - the previous style vector (before discovering v)
Thus the main loop below is for testing all offset positions in text sequencially.
The differences between currentStyle and v components (if found) are stored in style structure
as well as their positions, start and end.
So styles variable remembers all fragments positions with special formatting.
*/
function discoverStyles(text) {
  var styles = [[], [], [], [], [], []],
      currentStyle = [null, null, null, null, null, null],
      plainText = text.getText();

  for (var offset = 0; offset < plainText.length; offset++) {
    var v = [
      text.isBold(offset), text.isItalic(offset), text.isUnderline(offset),
      text.getForegroundColor(offset), text.getBackgroundColor(offset), text.getLinkUrl(offset)
    ];
    v.forEach(function(value, index) {
      if (value != currentStyle[index]) {  // Having a difference, we should..
        if (currentStyle[index] != null) { // close formatted fragment, which has been open before..
          styles[index][styles[index].length - 1].end = offset;
        }
        if (value != null) {  // then open new formatted fragment
          styles[index].push({start: offset, value: value});
        }
        currentStyle[index] = value;
      }
    });
  }
  return styles;
}
