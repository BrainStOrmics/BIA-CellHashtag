# Tool calling

[Frontend](/oss/python/langchain/frontend/overview)[Patterns](/oss/python/langchain/frontend/markdown-messages)

# Tool calling
Copy page

Display agent tool calls with rich, type-safe UI cardsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Agents can invoke external tools like weather APIs, calculators, web search,
database queries, and more. The results are in raw JSON. This pattern shows you
how to render
structured, type-safe UI cards for every tool call your agent makes, complete
with loading states and error handling.

## [​](#how-tool-calling-works)How tool calling works

When a LangGraph agent decides it needs external data, it emits one or more
**tool calls** as part of an AI message. Each tool call includes:

{' ' * (self.list_depth - 1)}- **name**: the tool being invoked (e.g. `"get_weather"`, `"calculator"`)

{' ' * (self.list_depth - 1)}- **args**: the structured arguments passed to the tool

{' ' * (self.list_depth - 1)}- **id**: a unique identifier linking the call to its result

The agent runtime executes the tool, and the result comes back as a
`ToolMessage`. The `useStream` hook unifies all of this into a single
`toolCalls` array you can render directly.

## [​](#setting-up-usestream)Setting up useStream

The first step is wiring up `useStream` to your agent backend. The hook returns
reactive state including a `toolCalls` array that updates in real time as the
agent streams.
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

const AGENT_URL = "http://localhost:2024";

export function Chat() {
 const stream = useStream<typeof myAgent>({
 apiUrl: AGENT_URL,
 assistantId: "tool_calling",
 });

 return (
 <div>
 {stream.messages.map((msg) => (
 <Message key={msg.id} message={msg} toolCalls={stream.toolCalls} />
 ))}
 </div>
 );
}

```

## [​](#the-toolcallwithresult-type)The ToolCallWithResult type

Each entry in the `toolCalls` array is a `ToolCallWithResult` object:

```
interface ToolCallWithResult {
 call: {
 id: string;
 name: string;
 args: Record<string, unknown>;
 };
 result: ToolMessage | undefined;
 state: "pending" | "completed" | "error";
}

```

| 
 | Property | Description
 | `call.id` | Unique ID matching the AI message’s `tool_calls` entry
 | `call.name` | The name of the tool (e.g. `"get_weather"`)
 | `call.args` | Structured arguments the agent passed to the tool
 | `result` | The `ToolMessage` response, available once the tool finishes
 | `state` | Lifecycle state: `"pending"` while running, `"completed"` on success, `"error"` on failure

## [​](#filtering-tool-calls-per-message)Filtering tool calls per message

An AI message may trigger multiple tool calls, and your chat may contain many AI
messages. To render the right tool cards under each message, filter by matching
`call.id` against the message’s `tool_calls` array:

```
function Message({
 message,
 toolCalls,
}: {
 message: AIMessage;
 toolCalls: ToolCallWithResult[];
}) {
 const messageToolCalls = toolCalls.filter((tc) =>
 message.tool_calls?.find((t) => t.id === tc.call.id)
 );

 return (
 <div>
 <p>{message.content}</p>
 {messageToolCalls.map((tc) => (
 <ToolCard key={tc.call.id} toolCall={tc} />
 ))}
 </div>
 );
}

```

## [​](#building-specialized-tool-cards)Building specialized tool cards

Rather than dumping raw JSON, build dedicated UI components for each tool. Use
`call.name` to select the right card:

```
function ToolCard({ toolCall }: { toolCall: ToolCallWithResult }) {
 if (toolCall.state === "pending") {
 return <LoadingCard name={toolCall.call.name} />;
 }

 if (toolCall.state === "error") {
 return <ErrorCard name={toolCall.call.name} error={toolCall.result} />;
 }

 switch (toolCall.call.name) {
 case "get_weather":
 return <WeatherCard args={toolCall.call.args} result={toolCall.result} />;
 case "calculator":
 return (
 <CalculatorCard args={toolCall.call.args} result={toolCall.result} />
 );
 case "web_search":
 return <SearchCard args={toolCall.call.args} result={toolCall.result} />;
 default:
 return <GenericToolCard toolCall={toolCall} />;
 }
}

```

### [​](#weather-card-example)Weather card example

```
function WeatherCard({
 args,
 result,
}: {
 args: { location: string };
 result: ToolMessage;
}) {
 const data = JSON.parse(result.content as string);

 return (
 <div className="rounded-lg border p-4">
 <div className="flex items-center gap-2">
 <CloudIcon />
 <h3 className="font-semibold">{args.location}</h3>
 </div>
 <div className="mt-2 text-3xl font-bold">{data.temperature}°F</div>
 <p className="text-muted-foreground">{data.condition}</p>
 </div>
 );
}

```

### [​](#loading-and-error-states)Loading and error states

Always handle the pending and error states to give users clear feedback:

```
function LoadingCard({ name }: { name: string }) {
 return (
 <div className="flex items-center gap-2 rounded-lg border p-4 animate-pulse">
 <Spinner />
 <span>Running {name}...</span>
 </div>
 );
}

function ErrorCard({ name, error }: { name: string; error?: ToolMessage }) {
 return (
 <div className="rounded-lg border border-red-300 bg-red-50 p-4">
 <h3 className="font-semibold text-red-700">Error in {name}</h3>
 <p className="text-sm text-red-600">
 {error?.content ?? "Tool execution failed"}
 </p>
 </div>
 );
}

```

## [​](#type-safe-tool-arguments)Type-safe tool arguments

If your tools are defined with structured schemas, you can use the
`ToolCallFromTool` utility type to get fully typed `args`:

```
import { tool } from "@langchain/core/tools";
import { z } from "zod";

const getWeather = tool(async ({ location }) => { /* ... */ }, {
 name: "get_weather",
 description: "Get the current weather for a location",
 schema: z.object({
 location: z.string().describe("City name"),
 }),
});

type WeatherToolCall = ToolCallFromTool<typeof getWeather>;
// WeatherToolCall.call.args is now { location: string }

```

Using `ToolCallFromTool` gives you compile-time safety. If the tool schema
changes, your UI components will flag type errors immediately.

## [​](#rendering-tool-calls-inline-with-streaming-text)Rendering tool calls inline with streaming text

Tool calls often arrive interleaved with streamed text. The `useStream` hook
keeps `toolCalls` in sync with the stream, so pending cards appear as soon as
the agent emits the call, before the tool has finished executing.
This means users see:

{' ' * (self.list_depth - 1)}- The AI’s text as it streams in

{' ' * (self.list_depth - 1)}- A loading card the moment a tool call is emitted

{' ' * (self.list_depth - 1)}- The card updates to show the result once the tool completes

Tool calls update in place. The same `call.id` transitions from `"pending"` to
`"completed"` (or `"error"`), so your UI re-renders the same component
with new state.

## [​](#handling-multiple-concurrent-tool-calls)Handling multiple concurrent tool calls

Agents can invoke several tools in parallel. The `toolCalls` array will contain
multiple entries with `state: "pending"` simultaneously. Each resolves
independently, so your UI should handle partial completion gracefully:

```
function ToolCallList({ toolCalls }: { toolCalls: ToolCallWithResult[] }) {
 const pending = toolCalls.filter((tc) => tc.state === "pending");
 const completed = toolCalls.filter((tc) => tc.state === "completed");

 return (
 <div className="space-y-2">
 {completed.map((tc) => (
 <ToolCard key={tc.call.id} toolCall={tc} />
 ))}
 {pending.map((tc) => (
 <LoadingCard key={tc.call.id} name={tc.call.name} />
 ))}
 </div>
 );
}

```

## [​](#best-practices)Best practices

Follow these guidelines when building tool call UIs:

{' ' * (self.list_depth - 1)}- **Always handle all three states**: `pending`, `completed`, and `error`.
Users should never see a blank card.

{' ' * (self.list_depth - 1)}- **Parse results safely**. Tool results arrive as strings. Wrap
`JSON.parse()` in a try/catch and show a fallback on parse failure.

{' ' * (self.list_depth - 1)}- **Provide a generic fallback**. Not every tool needs a bespoke card. Render
a collapsible JSON view for unknown tool names.

{' ' * (self.list_depth - 1)}- **Show the tool name and args during loading**. Users want to know *what*
the agent is doing, even before the result arrives.

{' ' * (self.list_depth - 1)}- **Keep cards compact**. Tool cards sit inline with chat messages. Avoid
overwhelming the conversation with oversized widgets.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/frontend/tool-calling.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Markdown messagesPrevious](/oss/python/langchain/frontend/markdown-messages)[Human-in-the-LoopNext](/oss/python/langchain/frontend/human-in-the-loop)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
