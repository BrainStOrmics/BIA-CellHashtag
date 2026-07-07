# Generative UI

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Generative UI
Copy page

Render AI-generated user interfaces using json-renderCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Generative UI lets the AI generate complete user interfaces from natural language
prompts. Instead of rendering text responses in chat bubbles, the AI output **is**
the UI: forms, cards, dashboards, and more. The developer defines which components
are available (the “catalog”), and the AI composes them into a valid UI tree.
This pattern uses [json-render](https://json-render.dev), the Generative UI framework,
to define component catalogs, generate specs with AI, and render them safely across
React, Vue, Svelte, and Angular.

## [​](#how-it-works)How it works

{' ' * (self.list_depth - 1)}- **Define a catalog**: declare what components the AI can use, with typed props

{' ' * (self.list_depth - 1)}- **Prompt the AI**: describe the UI you want in natural language

{' ' * (self.list_depth - 1)}- **AI generates a spec**: a JSON document describing the component tree

{' ' * (self.list_depth - 1)}- **Render safely**: json-render’s `Renderer` renders the spec using your components

The catalog acts as a guardrail: the AI can only use components you’ve defined,
with props that match your schema. The output is always predictable and safe.

## [​](#define-a-component-catalog)Define a component catalog

The catalog describes every component the AI is allowed to use. Each component has a
Zod schema for its props and a description that the AI reads to understand when to
use it:

```
import { defineCatalog } from "@json-render/core";
import { schema } from "@json-render/react/schema";
import { z } from "zod";

const catalog = defineCatalog(schema, {
 components: {
 Card: {
 description: "A card container with optional title and padding",
 props: z.object({
 title: z.string().optional(),
 padding: z.enum(["sm", "md", "lg"]).optional(),
 }),
 },
 TextInput: {
 description: "A text input field with optional label and placeholder",
 props: z.object({
 label: z.string().optional(),
 placeholder: z.string().optional(),
 type: z.enum(["text", "email", "password", "number", "textarea"]).optional(),
 }),
 },
 Button: {
 description: "A clickable button with label and style variants",
 props: z.object({
 label: z.string(),
 variant: z.enum(["primary", "secondary", "ghost", "link"]).optional(),
 fullWidth: z.boolean().optional(),
 }),
 },
 },
 actions: {},
});

```

Keep catalogs focused. Include only components the AI needs for the use case.
A smaller catalog produces better results than a kitchen-sink approach.

## [​](#build-a-component-registry)Build a component registry

The registry maps each catalog component to its actual rendering implementation.
Use `defineRegistry` to get type-safe bindings between the catalog props and
your component functions:
ReactVue
```
import { defineRegistry, Renderer, JSONUIProvider } from "@json-render/react";

const { registry } = defineRegistry(catalog, {
 components: {
 Card: ({ props, children }) => (
 <div className="card">
 {props.title && <h2>{props.title}</h2>}
 {children}
 </div>
 ),
 TextInput: ({ props }) => (
 <div>
 {props.label && <label>{props.label}</label>}
 <input type={props.type ?? "text"} placeholder={props.placeholder} />
 </div>
 ),
 Button: ({ props }) => (
 <button className={props.variant ?? "primary"}>
 {props.label}
 </button>
 ),
 },
});

```

## [​](#connect-to-the-agent)Connect to the agent

The agent uses structured output to return a json-render spec. Set up `useStream`
with your agent’s assistant ID, then extract the spec from the AI message’s
`tool_calls`:
ReactVueSvelteAngular
```
import { useStream } from "@langchain/react";
import { AIMessage } from "@langchain/core/messages";

function GenerativeUI() {
 const stream = useStream<typeof myAgent>({
 apiUrl: "http://localhost:2024",
 assistantId: "generative_ui",
 });

 const aiMessage = stream.messages.find(AIMessage.isInstance);
 const rawSpec = aiMessage?.tool_calls?.[0]?.args;

 // ... filter and render (see streaming section below)
}

```

## [​](#stream-and-render-progressively)Stream and render progressively

During streaming, the spec is built up incrementally. Elements arrive one at a
time and may initially lack `type` or `props`. Filter to only complete elements
and pass `loading={true}` to the `Renderer`, which tells it to silently skip
children that haven’t arrived yet. The UI builds up component by component:

```
/*
 * Filter the streamed spec to only include elements with valid type/props,
 * enabling progressive rendering as the AI response builds up. Passing
 * loading={true} to the Renderer tells it to skip missing children silently.
 */
const spec = (() => {
 if (!rawSpec?.root || !rawSpec?.elements) return null;
 const rootEl = rawSpec.elements[rawSpec.root];
 if (!rootEl?.type || rootEl?.props == null) return null;

 const safeElements = {};
 for (const [key, el] of Object.entries(rawSpec.elements)) {
 if (el?.type && el?.props != null) {
 safeElements[key] = el;
 }
 }
 return { root: rawSpec.root, elements: safeElements };
})();

return (
 <>
 {spec && (
 <JSONUIProvider registry={registry}>
 <Renderer spec={spec} registry={registry} loading={stream.isLoading} />
 </JSONUIProvider>
 )}
 </>
);

```

The `JSONUIProvider` is required to set up json-render’s internal context
providers (state, visibility, validation, actions). The `Renderer` component
must be rendered inside it.

## [​](#the-spec-format)The spec format

The AI agent generates a flat JSON spec with a `root` key pointing to the
root element and an `elements` map containing all components:

```
{
 "root": "login-card",
 "elements": {
 "login-card": {
 "type": "Card",
 "props": { "title": "Login" },
 "children": ["login-stack"]
 },
 "login-stack": {
 "type": "Stack",
 "props": { "direction": "vertical", "gap": "md" },
 "children": ["email-input", "password-input", "submit-btn"]
 },
 "email-input": {
 "type": "TextInput",
 "props": { "label": "Email", "placeholder": "Enter your email", "type": "email" },
 "children": []
 },
 "password-input": {
 "type": "TextInput",
 "props": { "label": "Password", "placeholder": "Enter your password", "type": "password" },
 "children": []
 },
 "submit-btn": {
 "type": "Button",
 "props": { "label": "Sign In", "variant": "primary", "fullWidth": true },
 "children": []
 }
 }
}

```

Each element references its children by ID, and leaf elements like `TextInput`
and `Button` have empty `children` arrays.

## [​](#best-practices)Best practices

{' ' * (self.list_depth - 1)}- **Use descriptive component descriptions**: the AI uses these to understand when
to use each component. Clear descriptions lead to better UI generation.

{' ' * (self.list_depth - 1)}- **Validate before rendering**: always check that elements have valid `type` and
non-null `props` before passing to the Renderer, since streaming delivers partial data.

{' ' * (self.list_depth - 1)}- **Design for streaming**: pass `loading={true}` during streaming so the Renderer
gracefully handles children that haven’t arrived yet. Users see the UI build up
in real time rather than waiting for the full response.

{' ' * (self.list_depth - 1)}- **Style with design tokens**: use CSS custom properties so rendered components
adapt to light and dark themes automatically.

{' ' * (self.list_depth - 1)}- **Wrap with JSONUIProvider**: the `Renderer` must be inside a `JSONUIProvider`
to access json-render’s internal context for state, visibility, and actions.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/generative-ui.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Time travelPrevious](/oss/python/langchain/frontend/time-travel)[OverviewNext](/oss/python/langchain/frontend/integrations/overview)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
