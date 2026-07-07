# Structured output

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Structured output
Copy page

Render structured agent responses with custom UI components instead of plain textCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Structured output lets the agent return typed, machine-readable data instead of plain text. Instead of rendering a single string, you get a structured object you can map to any UI: cards, tables, charts, step-by-step breakdowns, or domain-specific renderers.

## [​](#what-is-structured-output)What is structured output?

Instead of returning a free-form text response, the agent uses a tool call to return a structured object conforming to a predefined schema. This gives you:

{' ' * (self.list_depth - 1)}- **Type-safe data**: parse the response into a known TypeScript type

{' ' * (self.list_depth - 1)}- **Precise rendering control**: render each field with its own UI treatment

{' ' * (self.list_depth - 1)}- **Consistent formatting**: every response follows the same structure regardless of the underlying model

The agent accomplishes this by calling a “structured output” tool whose arguments contain the response data. The tool itself doesn’t execute any logic and is purely a vehicle for returning typed data.

## [​](#use-cases)Use cases

{' ' * (self.list_depth - 1)}- **Product comparisons**: feature tables, pros/cons lists, ratings

{' ' * (self.list_depth - 1)}- **Data analysis**: summaries with metrics, breakdowns, and highlights

{' ' * (self.list_depth - 1)}- **Step-by-step guides**: ordered instructions with descriptions and code snippets

{' ' * (self.list_depth - 1)}- **Recipes**: ingredients, steps, timings, and nutritional info

{' ' * (self.list_depth - 1)}- **Math and science**: formulas rendered with LaTeX, step-by-step derivations

{' ' * (self.list_depth - 1)}- **Travel planning**: itineraries with dates, locations, and cost estimates

## [​](#define-a-schema)Define a schema

Define a TypeScript type for the structured data the agent returns. The shape of this schema determines how you render the UI.
Here’s an example for a recipe assistant:

```
interface Ingredient {
 name: string;
 amount: string;
 unit: string;
}

interface RecipeStep {
 instruction: string;
 duration?: string;
}

interface Recipe {
 title: string;
 description: string;
 servings: number;
 ingredients: Ingredient[];
 steps: RecipeStep[];
 totalTime: string;
}

```

| 
 | Field | Type | Description
 | `title` | `string` | Name of the recipe
 | `description` | `string` | Short summary of the dish
 | `servings` | `number` | Number of servings
 | `ingredients` | `Ingredient[]` | List of ingredients with amounts and units
 | `steps` | `RecipeStep[]` | Ordered preparation steps
 | `totalTime` | `string` | Estimated total preparation and cooking time
Your schema can be anything. The pattern works the same way regardless of shape.

## [​](#extract-structured-output-from-messages)Extract structured output from messages

The structured output lives in the `tool_calls` array of the last `AIMessage`. Extract it by finding the AI message and accessing the first tool call’s arguments:

```
import { AIMessage } from "@langchain/core/messages";

function extractStructuredOutput<T>(messages: any[]): T | null {
 const aiMessages = messages.filter(AIMessage.isInstance);
 if (aiMessages.length === 0) return null;

 const lastAI = aiMessages[aiMessages.length - 1];
 const toolCall = lastAI.tool_calls?.[0];
 if (!toolCall) return null;

 return toolCall.args as T;
}

```

The structured output tool call may not have `args` populated until the agent finishes streaming. During streaming, `args` may be partially populated or undefined. Always check for completeness before rendering.

## [​](#set-up-usestream)Set up `useStream`

Define a TypeScript interface matching your agent’s state schema and pass it as a type parameter to `useStream` for type-safe access to state values. In the examples below, replace `typeof myAgent` with your interface name:

```
import type { BaseMessage } from "@langchain/core/messages";

interface AgentState {
 messages: BaseMessage[];
}

```

ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";
import { AIMessage } from "@langchain/core/messages";

function RecipeChat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "recipe_assistant",
 });

 const recipe = extractStructuredOutput<Recipe>(stream.messages);

 return (
 <div>
 {!recipe && !stream.isLoading && (
 <PromptInput onSubmit={(text) =>
 stream.submit({ messages: [{ type: "human", content: text }] })
 } />
 )}
 {stream.isLoading && <LoadingIndicator />}
 {recipe && <RecipeCard recipe={recipe} />}
 </div>
 );
}

```

## [​](#render-the-structured-data)Render the structured data

Once you have a typed object, build a component that maps each field to the
appropriate UI element. This is the core of the pattern: turning structured
data into a purpose-built interface.

```
function RecipeCard({ recipe }: { recipe: Recipe }) {
 return (
 <div className="recipe-card">
 <div className="recipe-header">
 <h3>{recipe.title}</h3>
 <p className="recipe-description">{recipe.description}</p>
 <div className="recipe-meta">
 <span>{recipe.servings} servings</span>
 <span>{recipe.totalTime}</span>
 </div>
 </div>

 <div className="recipe-ingredients">
 <h4>Ingredients</h4>
 <ul>
 {recipe.ingredients.map((ing, i) => (
 <li key={i}>
 <strong>{ing.amount} {ing.unit}</strong> {ing.name}
 </li>
 ))}
 </ul>
 </div>

 <div className="recipe-steps">
 <h4>Instructions</h4>
 {recipe.steps.map((step, i) => (
 <div key={i} className="step">
 <div className="step-number">Step {i + 1}</div>
 <p className="step-instruction">{step.instruction}</p>
 {step.duration && (
 <span className="step-duration">{step.duration}</span>
 )}
 </div>
 ))}
 </div>
 </div>
 );
}

```

The same approach works for any domain. Map each field to the UI element that best represents it:

| 
 | Data type | Rendering strategy
 | Plain text | Paragraphs, headings, list items
 | Numbers/metrics | Stat cards, progress bars, badges
 | Arrays | Lists, tables, grids
 | Nested objects | Nested cards, accordion sections
 | Markdown | Markdown renderer (e.g. `react-markdown`)
 | LaTeX/math | KaTeX or MathJax
 | Dates/times | Formatted timestamps, relative time
 | URLs | Links, embedded previews

## [​](#handle-partial-streaming-data)Handle partial streaming data

During streaming, the tool call arguments may be incomplete JSON. Guard against this in your extraction logic:

```
function extractStructuredOutput<T>(
 messages: any[],
 requiredFields: string[] = [],
): T | null {
 const aiMessages = messages.filter(AIMessage.isInstance);
 if (aiMessages.length === 0) return null;

 const lastAI = aiMessages[aiMessages.length - 1];
 const toolCall = lastAI.tool_calls?.[0];
 if (!toolCall?.args) return null;

 const args = toolCall.args as Record<string, unknown>;
 const hasRequired = requiredFields.every(
 (field) => args[field] !== undefined
 );

 if (requiredFields.length > 0 && !hasRequired) return null;
 return args as T;
}

```

Use the `requiredFields` parameter to wait until critical fields are populated before rendering:

```
const recipe = extractStructuredOutput<Recipe>(stream.messages, [
 "title",
 "ingredients",
 "steps",
]);

```

## [​](#render-progressively-during-streaming)Render progressively during streaming

Rather than waiting for the complete structured output, render fields as they arrive. This gives users immediate feedback while the agent is still generating:

```
function ProgressiveRecipeCard({ messages }: { messages: any[] }) {
 const partial = extractStructuredOutput<Partial<Recipe>>(messages);
 if (!partial) return null;

 return (
 <div className="recipe-card">
 {partial.title && <h3>{partial.title}</h3>}
 {partial.description && <p>{partial.description}</p>}

 {partial.ingredients && partial.ingredients.length > 0 && (
 <div className="recipe-ingredients">
 <h4>Ingredients</h4>
 <ul>
 {partial.ingredients.map((ing, i) => (
 <li key={i}>
 {ing.amount} {ing.unit} {ing.name}
 </li>
 ))}
 </ul>
 </div>
 )}

 {partial.steps && partial.steps.length > 0 && (
 <div className="recipe-steps">
 <h4>Instructions</h4>
 {partial.steps.map((step, i) => (
 <div key={i} className="step">
 <div className="step-number">Step {i + 1}</div>
 <p>{step.instruction}</p>
 </div>
 ))}
 </div>
 )}
 </div>
 );
}

```

Progressive rendering works well when the schema has a natural top-to-bottom
order: title, then description, then details. The agent typically generates
fields in schema order, so the UI fills in naturally.

## [​](#reset-and-re-submit)Reset and re-submit

To let the user submit a new query after viewing a result, add a button that starts a new thread:

```
{recipe && (
 <button onClick={() => stream.switchThread(null)}>
 Start over
 </button>
)}

```

This clears the current conversation and lets the user begin a fresh interaction.

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Validate before rendering**: always check that required fields exist before rendering, since streaming may deliver partial data

{' ' * (self.list_depth - 1)}- **Use a generic extraction function**: parameterize your extraction logic with a type and required fields so it works across different schemas

{' ' * (self.list_depth - 1)}- **Render progressively**: show fields as they arrive rather than waiting for the complete object, so users see immediate feedback

{' ' * (self.list_depth - 1)}- **Provide fallback representations**: if a field supports rich rendering (LaTeX, Markdown, charts), also include a plain-text equivalent in your schema as a fallback

{' ' * (self.list_depth - 1)}- **Keep schemas flat when possible**: deeply nested schemas are harder to render progressively and more likely to break during partial streaming

{' ' * (self.list_depth - 1)}- **Match UI to data**: choose the rendering strategy that best represents each field type (tables for arrays, cards for nested objects, badges for status fields)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/structured-output.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Reasoning tokensPrevious](/oss/python/langchain/frontend/reasoning-tokens)[Message queuesNext](/oss/python/langchain/frontend/message-queues)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
